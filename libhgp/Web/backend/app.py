"""
LibHGP Web 服务
后端：FastAPI
调用：hgp_py（pybind11 扩展模块）
hgp_py.cp312-win_amd64.pyd 位于 vis/ 根目录
启动方式：
cd vis/backend
python -m uvicorn app:app --reload --port 8000
"""

import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from fastapi import FastAPI, Response, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel
from collections import defaultdict
import hgp_py
import uuid

app = FastAPI(title="LibHGP Web Service")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

MESH_DIR   = os.path.join(os.path.dirname(__file__), "mesh_cache")
INDEX_HTML = os.path.join(os.path.dirname(__file__), "web", "index.html")
os.makedirs(MESH_DIR, exist_ok=True)

# ════════════════════════════════════════════════
# 静态页面
# ════════════════════════════════════════════════

@app.get("/")
async def root():
    return FileResponse(INDEX_HTML)

@app.get("/favicon.ico")
async def favicon():
    return Response(status_code=204)

# ════════════════════════════════════════════════
# Pydantic 请求模型（统一放这里，便于查阅）
# ════════════════════════════════════════════════

class PointsRequest(BaseModel):
    points: list[list[float]]

class MeshIdRequest(BaseModel):
    mesh_id: str

class SmoothRequest(BaseModel):
    vertices:   list[list[float]]
    face_id_0:  list[int]
    face_id_1:  list[int]
    face_id_2:  list[int]
    iterations: int = 5

class GeodesicRequest(BaseModel):
    mesh_id: str
    source:  list[float]
    target:  list[float]

class SlicerRequest(BaseModel):
    mesh_id:      str
    plane_normal: list[float]
    plane_ds:     list[float]

class CurvatureRequest(BaseModel):
    mesh_id:      str
    input_points: list[list[float]]

class SurfaceDecompRequest(BaseModel):
    mesh_id: str

class SamplingRequest(BaseModel):
    mesh_id: str
    density: float = 0.05

class PolygonRequest(BaseModel):
    polygon: list[list[float]]

class PolygonOffsetRequest(BaseModel):
    polygon:  list[list[float]]
    distance: float

class TwoPolygonRequest(BaseModel):
    polygon_0: list[list[float]]
    polygon_1: list[list[float]]

class Polygon2DSamplingRequest(BaseModel):
    polygon: list[list[float]]
    density: float = 0.05

class DistanceQueryRequest(BaseModel):
    mesh_id:      str
    query_points: list[list[float]]

class VisRequest(BaseModel):
    obj_string: str
    label:      str = ""

# ── CSG ──────────────────────────────────────────
class CSGHoleRequest(BaseModel):
    """在 mesh_id 对应的网格上，沿指定法线方向打一个圆柱孔"""
    mesh_id: str
    center:  list[float]   # 孔的圆心（网格表面上的拾取点）
    normal:  list[float]   # 孔的轴线方向（拾取面的法线）
    radius:  float
    depth:   float

class CSGBooleanRequest(BaseModel):
    """两个已上传的网格做通用布尔运算"""
    mesh_id_a: str
    mesh_id_b: str
    operation: str         # "difference" | "union" | "intersection"

class CSGPreviewRequest(BaseModel):
    """仅返回圆柱工具体网格，不执行布尔运算"""
    center: list[float]
    normal: list[float]
    radius: float
    depth:  float

# ════════════════════════════════════════════════
# 工具函数
# ════════════════════════════════════════════════

def mesh_response(vertices, fi0, fi1, fi2, mesh_id: str = None):
    result = {
        "vertices": [c for v in vertices for c in v],
        "faces":    [i for t in zip(fi0, fi1, fi2) for i in t],
    }
    if mesh_id:
        result["mesh_id"] = mesh_id
    return result


def repair_mesh(vertices, fi0, fi1, fi2):
    """
    修复网格，满足 CGAL 严格流形条件：
    1. 合并重复顶点  2. 删除退化面片  3. 删除重复面片
    4. 统一半边方向  5. 删除非流形边  6. 消除蝴蝶结顶点
    """
    EPS = 1e-8

    # ── 步骤 1：合并重复顶点 ─────────────────────
    unique_verts, old_to_new = [], {}
    for old_idx, v in enumerate(vertices):
        found = False
        for new_idx, uv in enumerate(unique_verts):
            if abs(v[0]-uv[0]) < EPS and abs(v[1]-uv[1]) < EPS and abs(v[2]-uv[2]) < EPS:
                old_to_new[old_idx] = new_idx
                found = True
                break
        if not found:
            old_to_new[old_idx] = len(unique_verts)
            unique_verts.append(v)

    # ── 步骤 2：删除退化面片 ─────────────────────
    new_fi0, new_fi1, new_fi2 = [], [], []
    for a, b, c in zip(fi0, fi1, fi2):
        na, nb, nc = old_to_new[a], old_to_new[b], old_to_new[c]
        if na == nb or nb == nc or na == nc:
            continue
        new_fi0.append(na); new_fi1.append(nb); new_fi2.append(nc)

    # ── 步骤 3：删除重复面片 ─────────────────────
    seen = set()
    dedup_fi0, dedup_fi1, dedup_fi2 = [], [], []
    for a, b, c in zip(new_fi0, new_fi1, new_fi2):
        key = tuple(sorted([a, b, c]))
        if key not in seen:
            seen.add(key)
            dedup_fi0.append(a); dedup_fi1.append(b); dedup_fi2.append(c)

    # ── 步骤 4：统一半边方向 ─────────────────────
    half_edge_map = {}
    fixed_fi0, fixed_fi1, fixed_fi2 = [], [], []
    for idx, (a, b, c) in enumerate(zip(dedup_fi0, dedup_fi1, dedup_fi2)):
        edges = [(a,b),(b,c),(c,a)]
        if any(e in half_edge_map for e in edges):
            a, b, c = a, c, b
            edges = [(a,b),(b,c),(c,a)]
        if any(e in half_edge_map for e in edges):
            continue
        for e in edges:
            half_edge_map[e] = idx
        fixed_fi0.append(a); fixed_fi1.append(b); fixed_fi2.append(c)

    # ── 步骤 5：删除非流形边 ─────────────────────
    edge_to_faces = defaultdict(list)
    for idx, (a, b, c) in enumerate(zip(fixed_fi0, fixed_fi1, fixed_fi2)):
        edge_to_faces[tuple(sorted([a,b]))].append(idx)
        edge_to_faces[tuple(sorted([b,c]))].append(idx)
        edge_to_faces[tuple(sorted([a,c]))].append(idx)

    bad_faces = set()
    for edge, faces in edge_to_faces.items():
        if len(faces) > 2:
            for f in faces[2:]: bad_faces.add(f)

    manifold_fi0, manifold_fi1, manifold_fi2 = [], [], []
    for idx, (a, b, c) in enumerate(zip(fixed_fi0, fixed_fi1, fixed_fi2)):
        if idx not in bad_faces:
            manifold_fi0.append(a); manifold_fi1.append(b); manifold_fi2.append(c)

    # ── 步骤 6：消除蝴蝶结顶点 ──────────────────
    def get_fan_components(v, face_indices, mfi0, mfi1, mfi2):
        if not face_indices: return []
        nbv_to_f = defaultdict(list)
        for f in face_indices:
            for u in (mfi0[f], mfi1[f], mfi2[f]):
                if u != v: nbv_to_f[u].append(f)
        adj = defaultdict(set)
        for u, fs in nbv_to_f.items():
            if len(fs) == 2:
                adj[fs[0]].add(fs[1]); adj[fs[1]].add(fs[0])
        visited, components = set(), []
        for start in face_indices:
            if start in visited: continue
            comp, queue = set(), [start]
            while queue:
                cur = queue.pop()
                if cur in visited: continue
                visited.add(cur); comp.add(cur)
                for nb in adj[cur]:
                    if nb not in visited: queue.append(nb)
            components.append(comp)
        return components

    v2f = defaultdict(list)
    for idx, (a, b, c) in enumerate(zip(manifold_fi0, manifold_fi1, manifold_fi2)):
        v2f[a].append(idx); v2f[b].append(idx); v2f[c].append(idx)

    bowtie_bad = set()
    for v, face_indices in v2f.items():
        comps = get_fan_components(v, face_indices, manifold_fi0, manifold_fi1, manifold_fi2)
        if len(comps) > 1:
            largest = max(comps, key=len)
            for comp in comps:
                if comp is not largest: bowtie_bad |= comp

    final_fi0, final_fi1, final_fi2 = [], [], []
    for idx, (a, b, c) in enumerate(zip(manifold_fi0, manifold_fi1, manifold_fi2)):
        if idx not in bowtie_bad:
            final_fi0.append(a); final_fi1.append(b); final_fi2.append(c)

    print(f"[repair_mesh] {len(vertices)}v {len(fi0)}f → "
          f"{len(unique_verts)}v {len(final_fi0)}f "
          f"(非流形={len(bad_faces)}, 蝴蝶结={len(bowtie_bad)})")
    return unique_verts, final_fi0, final_fi1, final_fi2

def _normalize_mesh(vertices, fi0, fi1, fi2,
                    target_size: float = 2.0):
    """
    将网格平移到原点并等比缩放至目标尺寸。

    target_size : 归一化后包围盒最长边的长度（默认 2.0，
                  对应 Three.js 场景里约 [-1, 1] 的舒适视野）

    返回变换后的 (vertices, fi0, fi1, fi2)（fi* 不变）
    """
    if not vertices:
        return vertices, fi0, fi1, fi2

    xs = [v[0] for v in vertices]
    ys = [v[1] for v in vertices]
    zs = [v[2] for v in vertices]

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    # 中心点
    cx = (min_x + max_x) / 2.0
    cy = (min_y + max_y) / 2.0
    cz = (min_z + max_z) / 2.0

    # 等比缩放系数（最长边 → target_size）
    span = max(max_x - min_x, max_y - min_y, max_z - min_z)
    scale = (target_size / span) if span > 1e-12 else 1.0

    new_verts = [
        [(v[0] - cx) * scale,
         (v[1] - cy) * scale,
         (v[2] - cz) * scale]
        for v in vertices
    ]
    return new_verts, fi0, fi1, fi2

def _center_mesh(vertices, fi0, fi1, fi2):
    """
    将网格包围盒中心平移到原点。
    仅做平移，不缩放——保留所有实际距离/角度信息。
    返回平移后的 (vertices, fi0, fi1, fi2)（fi* 不变）。
    """
    if not vertices:
        return vertices, fi0, fi1, fi2

    xs = [v[0] for v in vertices]
    ys = [v[1] for v in vertices]
    zs = [v[2] for v in vertices]

    cx = (min(xs) + max(xs)) / 2.0
    cy = (min(ys) + max(ys)) / 2.0
    cz = (min(zs) + max(zs)) / 2.0

    new_verts = [[v[0] - cx, v[1] - cy, v[2] - cz] for v in vertices]
    return new_verts, fi0, fi1, fi2

def save_mesh(vertices, fi0, fi1, fi2, repair: bool = False) -> tuple[str, list, list, list, list]:
    """返回 (mesh_id, 最终顶点, fi0, fi1, fi2)，保证与磁盘文件一致"""
    if repair:
        vertices, fi0, fi1, fi2 = repair_mesh(vertices, fi0, fi1, fi2)
    mesh_id  = uuid.uuid4().hex[:8]
    obj_path = os.path.join(MESH_DIR, f"{mesh_id}.obj")
    off_path = os.path.join(MESH_DIR, f"{mesh_id}.off")
    with open(obj_path, "w") as f:
        f.write("# LibHGP mesh cache\n")
        for v in vertices:
            f.write(f"v {v[0]} {v[1]} {v[2]}\n")
        for i in range(len(fi0)):
            f.write(f"f {fi0[i]+1} {fi1[i]+1} {fi2[i]+1}\n")
    nv, nf = len(vertices), len(fi0)
    with open(off_path, "w") as f:
        f.write("OFF\n")
        f.write(f"{nv} {nf} 0\n")
        for v in vertices:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")
        for i in range(nf):
            f.write(f"3 {fi0[i]} {fi1[i]} {fi2[i]}\n")
    return mesh_id, vertices, fi0, fi1, fi2


def get_mesh_path(mesh_id: str, ext: str = "obj") -> str:
    path = os.path.join(MESH_DIR, f"{mesh_id}.{ext}")
    if not os.path.exists(path):
        raise HTTPException(status_code=404,
            detail=f"mesh_id={mesh_id} 对应的 .{ext} 文件不存在")
    return path


def _load_mesh_from_obj(path: str):
    """读取 OBJ 文件，返回 (verts, fi0, fi1, fi2)"""
    verts, fi0, fi1, fi2 = [], [], [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts: continue
            if parts[0] == 'v':
                verts.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif parts[0] == 'f':
                fi0.append(int(parts[1]) - 1)
                fi1.append(int(parts[2]) - 1)
                fi2.append(int(parts[3]) - 1)
    return verts, fi0, fi1, fi2


def parse_obj_string(obj_str: str):
    """从 OBJ 文本字符串解析，返回 (verts, fi0, fi1, fi2)"""
    verts, fi0, fi1, fi2 = [], [], [], []
    def parse_idx(token): return int(token.split('/')[0]) - 1
    for line in obj_str.splitlines():
        parts = line.split()
        if not parts: continue
        if parts[0] == 'v' and len(parts) >= 4:
            verts.append([float(parts[1]), float(parts[2]), float(parts[3])])
        elif parts[0] == 'f' and len(parts) >= 4:
            fi0.append(parse_idx(parts[1]))
            fi1.append(parse_idx(parts[2]))
            fi2.append(parse_idx(parts[3]))
    return verts, fi0, fi1, fi2


# ── CSG 工具：归一化法线 ────────────────────────

def _normalize(v: list[float]) -> list[float]:
    nx, ny, nz = v
    length = math.sqrt(nx*nx + ny*ny + nz*nz)
    if length < 1e-9:
        raise HTTPException(status_code=400, detail="法线向量长度为零")
    return [nx/length, ny/length, nz/length]


def _do_boolean(verts_a, fi0_a, fi1_a, fi2_a,
                verts_b, fi0_b, fi1_b, fi2_b,
                operation: str):
    """调用 hgp_py CSG 布尔运算，返回 (verts, fi0, fi1, fi2)"""
    if operation == "union":
        ok, verts_out, fi0_out, fi1_out, fi2_out = hgp_py.HGP_Mesh_CSG_Union(
            verts_a, fi0_a, fi1_a, fi2_a,
            verts_b, fi0_b, fi1_b, fi2_b)
    elif operation == "difference":
        ok, verts_out, fi0_out, fi1_out, fi2_out = hgp_py.HGP_Mesh_CSG_Difference(
            verts_a, fi0_a, fi1_a, fi2_a,
            verts_b, fi0_b, fi1_b, fi2_b)
    elif operation == "intersection":
        ok, verts_out, fi0_out, fi1_out, fi2_out = hgp_py.HGP_Mesh_CSG_Intersection(
            verts_a, fi0_a, fi1_a, fi2_a,
            verts_b, fi0_b, fi1_b, fi2_b)
    else:
        raise HTTPException(status_code=400, detail=f"未知操作: {operation}，"
                            "合法值为 difference / union / intersection")
    if not ok or not verts_out:
        raise HTTPException(status_code=422,
            detail="布尔运算失败或结果为空，请确认两个网格均为封闭流形且存在交叠区域。")
    return verts_out, fi0_out, fi1_out, fi2_out

# ════════════════════════════════════════════════
# 接口 1：3D 凸包
# ════════════════════════════════════════════════

@app.post("/api/convex_hull_3d")
def convex_hull_3d(req: PointsRequest):
    hull_verts, fi0, fi1, fi2 = hgp_py.HGP_3D_Convex_Hulls_C2(req.points)
    mesh_id, hull_verts, fi0, fi1, fi2 = save_mesh(hull_verts, fi0, fi1, fi2, repair=False)
    return mesh_response(hull_verts, fi0, fi1, fi2, mesh_id)

# ════════════════════════════════════════════════
# 接口 2：上传 OBJ 文件
# ════════════════════════════════════════════════

@app.post("/api/upload_obj")
async def upload_obj(file: UploadFile = File(...)):
    content = await file.read()
    text = content.decode("utf-8", errors="ignore")
    verts, fi0, fi1, fi2 = parse_obj_string(text)
    if not verts or not fi0:
        raise HTTPException(status_code=400, detail="OBJ 文件解析失败或不含三角面片")
    #verts, fi0, fi1, fi2 = _normalize_mesh(verts, fi0, fi1, fi2, target_size=2.0) # 归一化到 [-1,1] 范围，适合直接预览，但会丢失原始尺寸信息
    verts, fi0, fi1, fi2 = _center_mesh(verts, fi0, fi1, fi2)  # 仅平移到原点，保留尺寸信息，适合后续计算（如测地线）使用
    mesh_id, verts, fi0, fi1, fi2 = save_mesh(verts, fi0, fi1, fi2, repair=True)
    return mesh_response(verts, fi0, fi1, fi2, mesh_id)

# ════════════════════════════════════════════════
# 接口 3：OBJ 字符串直接渲染
# ════════════════════════════════════════════════

@app.post("/api/visualize")
def visualize(req: VisRequest):
    verts, fi0, fi1, fi2 = parse_obj_string(req.obj_string)
    mesh_id, verts, fi0, fi1, fi2 = save_mesh(verts, fi0, fi1, fi2, repair=False)
    return mesh_response(verts, fi0, fi1, fi2, mesh_id)

# ════════════════════════════════════════════════
# 接口 4：拉普拉斯平滑
# ════════════════════════════════════════════════

@app.post("/api/smooth")
def smooth(req: SmoothRequest):
    verts_out, fi0, fi1, fi2 = hgp_py.HGP_Mesh_Laplace_Smooth_C2(
        req.vertices, req.face_id_0, req.face_id_1, req.face_id_2, req.iterations)
    mesh_id, verts_out, fi0, fi1, fi2 = save_mesh(verts_out, fi0, fi1, fi2, repair=False)
    return mesh_response(verts_out, fi0, fi1, fi2, mesh_id)

# ════════════════════════════════════════════════
# 接口 5：Loop 细分
# ════════════════════════════════════════════════

@app.post("/api/subdivide")
def subdivide(req: SmoothRequest):
    verts = [list(v) for v in req.vertices]
    fi0   = list(req.face_id_0)
    fi1   = list(req.face_id_1)
    fi2   = list(req.face_id_2)
    hgp_py.HGP_Mesh_Loop_Subdivision_One_Step(verts, fi0, fi1, fi2)
    mesh_id, verts, fi0, fi1, fi2 = save_mesh(verts, fi0, fi1, fi2, repair=False)
    return mesh_response(verts, fi0, fi1, fi2, mesh_id)

# ════════════════════════════════════════════════
# 接口 6：曲率热力图
# ════════════════════════════════════════════════

@app.post("/api/curvature")
def curvature(req: CurvatureRequest):
    mesh_path = get_mesh_path(req.mesh_id)
    max_curs, min_curs, max_dirs, min_dirs = hgp_py.HGP_Curvature_Mesh(
        mesh_path, req.input_points)
    return {
        "max_curvatures": max_curs,
        "min_curvatures": min_curs,
        "max_directions": [[d[0],d[1],d[2]] for d in max_dirs],
        "min_directions": [[d[0],d[1],d[2]] for d in min_dirs],
    }

# ════════════════════════════════════════════════
# 接口 7：最短测地路径 & 距离
# ════════════════════════════════════════════════

@app.post("/api/geodesic_path")
def geodesic_path(req: GeodesicRequest):
    mesh_path   = get_mesh_path(req.mesh_id)
    path_points = hgp_py.HGP_Shortest_Geodesic_Path_C3(mesh_path, req.source, req.target)
    if len(path_points) < 2:
        raise HTTPException(status_code=422,
            detail="两点不连通：source 和 target 必须在同一连通网格上。")
    dist = hgp_py.HGP_Geodesic_Distance(mesh_path, req.source, req.target)
    if dist < 0:
        raise HTTPException(status_code=422, detail="测地线计算失败。")
    return {"path": [[p[0],p[1],p[2]] for p in path_points], "distance": dist}

@app.post("/api/geodesic_distance")
def geodesic_distance(req: GeodesicRequest):
    mesh_path = get_mesh_path(req.mesh_id)
    dist = hgp_py.HGP_Geodesic_Distance(mesh_path, req.source, req.target)
    return {"distance": dist}

# ════════════════════════════════════════════════
# 接口 8：网格切片
# ════════════════════════════════════════════════

@app.post("/api/slicer")
def slicer(req: SlicerRequest):
    mesh_path = get_mesh_path(req.mesh_id, ext="off")
    offsetses, offsets = hgp_py.HGP_Slicer_Mesh(
        mesh_path, req.plane_normal, req.plane_ds
    )
    print(f"plane_ds: {req.plane_ds}")
    print(f"offsetses 长度: {len(offsetses)}, 各平面段数: {[len(p) for p in offsetses]}")
    print(f"offsets 长度: {len(offsets)}")
    return {
        "slices": [
            {
                "contours": [
                    [[p[0], p[1], p[2]] for p in seg]
                    for seg in plane_contours
                ]
            }
            for plane_contours in offsetses
        ]
    }

# ════════════════════════════════════════════════
# 接口 9：网格边界提取
# ══════════════════════════════════════��═════════

@app.post("/api/mesh_boundary")
def mesh_boundary(req: MeshIdRequest):
    mesh_path  = get_mesh_path(req.mesh_id, ext="off")
    boundaries = hgp_py.HGP_3D_Triangle_Mesh_Boundary_C5(mesh_path)
    return {"boundaries": [[[p[0],p[1],p[2]] for p in loop] for loop in boundaries]}

# ════════════════════════════════════════════════
# 接口 10：顶点法向量
# ════════════════════════════════════════════════

@app.post("/api/mesh_normals")
def mesh_normals(req: CurvatureRequest):
    mesh_path = get_mesh_path(req.mesh_id)
    verts, fi0, fi1, fi2 = _load_mesh_from_obj(mesh_path)
    normals = hgp_py.HGP_3D_Mesh_Normal_C2(verts, fi0, fi1, fi2)
    return {
        "vertices": [[v[0],v[1],v[2]] for v in verts],
        "normals":  [[n[0],n[1],n[2]] for n in normals],
    }

# ════════════════════════════════════════════════
# 接口 11：包围盒
# ════════════════════════════════════════════════

@app.post("/api/mesh_bbox")
def mesh_bbox(req: CurvatureRequest):
    mesh_path = get_mesh_path(req.mesh_id)
    verts, fi0, fi1, fi2 = _load_mesh_from_obj(mesh_path)
    min_c, max_c = hgp_py.HGP_3D_Mesh_Boundingbox_C2(verts)
    return {
        "min_corner": [min_c[0],min_c[1],min_c[2]],
        "max_corner": [max_c[0],max_c[1],max_c[2]],
    }

# ════════════════════════════════════════════════
# 接口 12：表面采样
# ════════════════════════════════════════════════

@app.post("/api/mesh_sampling")
def mesh_sampling(req: SamplingRequest):
    mesh_path = get_mesh_path(req.mesh_id)
    pts = hgp_py.HGP_3D_Mesh_Regular_Sampling_C1(mesh_path, req.density)
    return {"points": [[p[0],p[1],p[2]] for p in pts]}

# ════════════════════════════════════════════════
# 接口 13：表面分解
# ════════════════════════════════════════════════

@app.post("/api/surface_decomposition")
def surface_decomposition(req: SurfaceDecompRequest):
    mesh_path = get_mesh_path(req.mesh_id)
    face_sdf, regions_nb, face_segments = hgp_py.HGP_Surface_Decomposition(mesh_path)
    return {"regions_count": regions_nb, "face_sdf": face_sdf, "face_segments": face_segments}

# ════════════════════════════════════════════════
# 接口 14：最近点到网格距离
# ════════════════════════════════════════════════

@app.post("/api/point_mesh_distance")
def point_mesh_distance(req: DistanceQueryRequest):
    mesh_path = get_mesh_path(req.mesh_id)
    distances = hgp_py.HGP_3D_Distance_Point_Mesh(mesh_path, req.query_points)
    return {"distances": distances}

# ════════════════════════════════════════════════
# 接口 15：2D 凸包
# ════════════════════════════════════════════════

@app.post("/api/convex_hull_2d")
def convex_hull_2d(req: PolygonRequest):
    hull = hgp_py.HGP_2D_Convex_Hulls(req.polygon)
    return {"hull": hull}

# ════════════════════════════════════════════════
# 接口 16：2D 多边形 Offset
# ════════════════════════════════════════════════

@app.post("/api/polygon_offset_2d")
def polygon_offset_2d(req: PolygonOffsetRequest):
    result_polys = hgp_py.HGP_2D_Polygon_One_Offsets(req.polygon, req.distance)
    return {"polygons": result_polys}

# ════════════════════════════════════════════════
# 接口 17：2D 多边形布尔并集
# ════════════════════════════════════════════════

@app.post("/api/polygon_union_2d")
def polygon_union_2d(req: TwoPolygonRequest):
    area, union_polys = hgp_py.HGP_2D_Two_Polygons_Union(req.polygon_0, req.polygon_1)
    return {"area": area, "polygons": union_polys}

# ════════════════════════════════════════════════
# 接口 18：2D 多边形均匀采样
# ════════════════════════════════════════════════

@app.post("/api/polygon_sampling_2d")
def polygon_sampling_2d(req: PolygonRequest, density: float = 0.05):
    pts = hgp_py.HGP_2D_Polygon_Regular_Sampling_C1(req.polygon, density)
    return {"points": pts}

@app.post("/api/polygon_sampling_2d_regular")
def polygon_sampling_2d_regular(req: Polygon2DSamplingRequest):
    pts = hgp_py.HGP_2D_Polygon_Regular_Sampling_C1(req.polygon, req.density)
    return {"points": pts}

# ════════════════════════════════════════════════
# 接口 19-21：CSG（构造实体几何）
# ════════════════════════════════════════════════
def _make_oriented_cylinder(center, normal_normalized, radius,
                             inward_ext, outward_ext,
                             segments=48):
    """
    注意：+normal方向=outward，-normal方向=inward
    生成圆柱，以 center 为拾取点（入口面）：
      - 朝 -normal（入射侧，模型内部）延伸 inward_ext
      - 朝 +normal（出射侧，摄像头方向）延伸 outward_ext

    坐标推导：
      生成底面在原点（z=0）、顶面在 z=total_height 的圆柱
      → 绕原点旋转对齐 normal（底面原点不动，顶面跑到 +normal*total_height）
      → 平移 (center - normal*outward_ext)
      → 底面最终在 center - normal*outward_ext（朝 -normal 方向，即入射端，深入模型）
      → 顶面最终在 center + normal*inward_ext（朝 +normal 方向，即出射端，超出表面）

    注意：参数名 inward/outward 与实际方向相反（历史遗留）
    调用时须遵守如下约定：
      csg_hole:             outward_ext=shortest_side, inward_ext=depth
      csg_preview_cylinder: outward_ext=depth,         inward_ext=0
    """
    ax, ay, az = normal_normalized
    total_height = inward_ext + outward_ext

    # ── 1. 生成底面在原点、顶面在 z=total_height 的圆柱 ──
    # HGP_Mesh_Make_Cylinder(cx,cy,cz, r, h, seg) 中 cz 是轴向中心
    # 轴向中心在 z=total_height/2 → 传入 cz=total_height/2
    half_total = total_height / 2.0
    cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2 = hgp_py.HGP_Mesh_Make_Cylinder(
        0.0, 0.0, half_total,
        radius, total_height, segments)
    # 产生：底面在 half_total - half_total = 0，顶面在 half_total + half_total = total_height ✓

    # ── 2. 旋转对齐 normal（绕原点旋转，底面圆心在原点保持不动）──
    z_axis = [0.0, 0.0, 1.0]
    dot = max(-1.0, min(1.0, ax*z_axis[0] + ay*z_axis[1] + az*z_axis[2]))

    if abs(dot - 1.0) > 1e-6:
        if abs(dot + 1.0) < 1e-6:
            # normal 与 +Z 完全反向，叉积为零，绕 X 轴旋转 π
            rot_axis = [1.0, 0.0, 0.0]
            angle = math.pi
        else:
            rx = z_axis[1]*az - z_axis[2]*ay
            ry = z_axis[2]*ax - z_axis[0]*az
            rz = z_axis[0]*ay - z_axis[1]*ax
            rl = math.sqrt(rx*rx + ry*ry + rz*rz)
            rot_axis = [rx/rl, ry/rl, rz/rl]
            angle = math.acos(dot)   # 弧度制

        tmp_id, _, _, _, _ = save_mesh(cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2, repair=False)
        tmp_obj = get_mesh_path(tmp_id, ext="obj")
        hgp_py.HGP_Rotation_Obj(tmp_obj, angle, rot_axis)
        cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2 = _load_mesh_from_obj(tmp_obj)

    # ── 3. 平移：底面（原点）→ center - normal * outward_ext ──
    ox = center[0] - ax * outward_ext
    oy = center[1] - ay * outward_ext
    oz = center[2] - az * outward_ext
    cyl_verts = [[v[0]+ox, v[1]+oy, v[2]+oz] for v in cyl_verts]

    return cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2

def is_close_by_magnitude(depth, thickness, ratio=1e-1):
    """
    根据 depth 的数量级决定绝对误差阈值
    ratio: 相对当前数量级的比例，默认 1%
    """
    # 获取 depth 的数量级（10 的幂次）
    if depth == 0:
        magnitude = 1  # 避免 log10(0)
    else:
        magnitude = 10 ** int(math.log10(abs(depth)))
    
    # 绝对误差阈值 = 数量级 * ratio
    tolerance = magnitude * ratio
    
    return abs(depth - thickness) <= tolerance

@app.post("/api/csg/hole")
def csg_hole(req: CSGHoleRequest):
    """在主体网格上打一个圆柱孔（difference 布尔运算）"""
    ax, ay, az = _normalize(req.normal)
    mesh_path = get_mesh_path(req.mesh_id)
    off_path  = get_mesh_path(req.mesh_id, ext="off")  # .off，用于 Surface_mesh 射线检测
    # ── 包围盒对角线 ──────────────────────────────
    verts_a, fi0_a, fi1_a, fi2_a = _load_mesh_from_obj(mesh_path)
    min_c, max_c = hgp_py.HGP_3D_Mesh_Boundingbox_C2(verts_a)
    dx = max_c[0]-min_c[0]; dy = max_c[1]-min_c[1]; dz = max_c[2]-min_c[2]
    diag = math.sqrt(dx*dx + dy*dy + dz*dz)

    # ── 出射侧：无条件延伸至 diag ─────────────────
    outward_ext = diag*0.01

    # ── 入射侧：射线检测通孔/盲孔 ────────────────
    EPS_OFFSET = max(diag * 0.005, 1e-4)
    shifted = [req.center[0] - ax*EPS_OFFSET,
               req.center[1] - ay*EPS_OFFSET,
               req.center[2] - az*EPS_OFFSET]
    ray_dir = [-ax, -ay, -az]
    hits = hgp_py.HGP_3D_Intersection_Rays_Mesh_Vector3d(
        [shifted], [ray_dir], off_path)
    print(f"hits: {hits}")
    print(f"hits type: {type(hits)}")
    print(f"hits length: {len(hits) if hits else 'None/Empty'}")
    print(f"shifted: {shifted}")
    thickness = None
    if hits:
        h = hits[0]
        ddx = h[0]-shifted[0]; ddy = h[1]-shifted[1]; ddz = h[2]-shifted[2]
        dist = math.sqrt(ddx*ddx + ddy*ddy + ddz*ddz)
        print(f"dist: {dist}")
        if dist > EPS_OFFSET * 0.1:
            thickness = EPS_OFFSET + dist

    if thickness is not None and is_close_by_magnitude(req.depth, thickness):
        small_eps = max(diag * 0.02, thickness * 0.05)
        inward_ext = req.depth + small_eps
        print(f"[csg_hole] 通孔: depth={req.depth:.4f} 约等于 thickness={thickness:.4f}")
    else:
        inward_ext = req.depth
        if thickness is None:
            print(f"[csg_hole] 射线未命中对面，按盲孔处理")
        else:
            print(f"[csg_hole] 盲孔: depth={req.depth:.4f} < thickness={thickness:.4f}")

    cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2 = _make_oriented_cylinder(
        req.center, ray_dir, req.radius,
        inward_ext=inward_ext,
        outward_ext=outward_ext)

    verts_out, fi0_out, fi1_out, fi2_out = _do_boolean(
        verts_a, fi0_a, fi1_a, fi2_a,
        cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2,
        "difference")

    mesh_id, verts_out, fi0_out, fi1_out, fi2_out = save_mesh(
        verts_out, fi0_out, fi1_out, fi2_out, repair=False)
    return mesh_response(verts_out, fi0_out, fi1_out, fi2_out, mesh_id)


@app.post("/api/csg/boolean")
def csg_boolean(req: CSGBooleanRequest):
    """对两个已上传的网格执行通用布尔运算"""
    verts_a, fi0_a, fi1_a, fi2_a = _load_mesh_from_obj(get_mesh_path(req.mesh_id_a))
    verts_b, fi0_b, fi1_b, fi2_b = _load_mesh_from_obj(get_mesh_path(req.mesh_id_b))
    verts_out, fi0_out, fi1_out, fi2_out = _do_boolean(
        verts_a, fi0_a, fi1_a, fi2_a,
        verts_b, fi0_b, fi1_b, fi2_b,
        req.operation)
    mesh_id, verts_out, fi0_out, fi1_out, fi2_out = save_mesh(verts_out, fi0_out, fi1_out, fi2_out, repair=False)
    return mesh_response(verts_out, fi0_out, fi1_out, fi2_out, mesh_id)


@app.post("/api/csg/preview_cylinder")
def csg_preview_cylinder(req: CSGPreviewRequest):
    """返回工具圆柱网格供前端半透明预览
    预览圆柱与实际打孔圆柱完全一致：
      - inward_ext = depth（朝模型内钻入深度）
      - outward_ext = 0（出射侧不超出，预览只展示钻入部分）
    """
    ax, ay, az = _normalize(req.normal)
    ray_dir = [-ax, -ay, -az]
    cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2 = _make_oriented_cylinder(
        req.center, ray_dir, req.radius,
        inward_ext=req.depth,
        outward_ext=0.0)
    mesh_id, cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2 = save_mesh(
        cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2, repair=False)
    return mesh_response(cyl_verts, cyl_fi0, cyl_fi1, cyl_fi2, mesh_id)
# ════════════════════════════════════════════════
# 接口：按 mesh_id 重新加载网格（用于撤销）
# ════════════════════════════════════════════════

@app.post("/api/mesh/reload")
def mesh_reload(req: MeshIdRequest):
    """按 mesh_id 从缓存重新读取并返回网格数据（供前端撤销使用）"""
    mesh_path = get_mesh_path(req.mesh_id)
    verts, fi0, fi1, fi2 = _load_mesh_from_obj(mesh_path)
    return mesh_response(verts, fi0, fi1, fi2, req.mesh_id)