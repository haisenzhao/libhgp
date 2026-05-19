// hgp_py.cpp — 完整绑定 libhgp.h 所有函数
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "libhgp.h"

namespace py = pybind11;

// ============================================================
// 辅助类型构造函数
// ============================================================
static inline Vector2d  to_vec2d (const std::vector<double>& v) {
    if (v.size()<2) throw std::runtime_error("need >=2 doubles for Vector2d");
    return Vector2d(v[0],v[1]);
}
static inline Vector3d  to_vec3d (const std::vector<double>& v) {
    if (v.size()<3) throw std::runtime_error("need >=3 doubles for Vector3d");
    return Vector3d(v[0],v[1],v[2]);
}
static inline std::vector<double> from_vec2d(const Vector2d& v){ return {v.x,v.y}; }
static inline std::vector<double> from_vec3d(const Vector3d& v){ return {v.x,v.y,v.z}; }

static inline Vector2d1 to_vec2d1(const std::vector<std::vector<double>>& vs){
    Vector2d1 r; r.reserve(vs.size());
    for(auto& v:vs) r.push_back(to_vec2d(v)); return r;
}
static inline std::vector<std::vector<double>> from_vec2d1(const Vector2d1& vs){
    std::vector<std::vector<double>> r; r.reserve(vs.size());
    for(auto& v:vs) r.push_back(from_vec2d(v)); return r;
}
static inline Vector2d2 to_vec2d2(const std::vector<std::vector<std::vector<double>>>& vss){
    Vector2d2 r; r.reserve(vss.size());
    for(auto& vs:vss) r.push_back(to_vec2d1(vs)); return r;
}
static inline std::vector<std::vector<std::vector<double>>> from_vec2d2(const Vector2d2& vss){
    std::vector<std::vector<std::vector<double>>> r; r.reserve(vss.size());
    for(auto& vs:vss) r.push_back(from_vec2d1(vs)); return r;
}
static inline Vector3d1 to_vec3d1(const std::vector<std::vector<double>>& vs){
    Vector3d1 r; r.reserve(vs.size());
    for(auto& v:vs) r.push_back(to_vec3d(v)); return r;
}
static inline std::vector<std::vector<double>> from_vec3d1(const Vector3d1& vs){
    std::vector<std::vector<double>> r; r.reserve(vs.size());
    for(auto& v:vs) r.push_back(from_vec3d(v)); return r;
}
static inline Vector3d2 to_vec3d2(const std::vector<std::vector<std::vector<double>>>& vss){
    Vector3d2 r; r.reserve(vss.size());
    for(auto& vs:vss) r.push_back(to_vec3d1(vs)); return r;
}
static inline std::vector<std::vector<std::vector<double>>> from_vec3d2(const Vector3d2& vss){
    std::vector<std::vector<std::vector<double>>> r; r.reserve(vss.size());
    for(auto& vs:vss) r.push_back(from_vec3d1(vs)); return r;
}
static inline std::vector<std::vector<std::vector<std::vector<double>>>> from_vec3d3(const Vector3d3& vsss){
    std::vector<std::vector<std::vector<std::vector<double>>>> r; r.reserve(vsss.size());
    for(auto& vss:vsss) r.push_back(from_vec3d2(vss)); return r;
}
static inline Vector3d3 to_vec3d3(const std::vector<std::vector<std::vector<std::vector<double>>>>& vsss){
    Vector3d3 r; r.reserve(vsss.size());
    for(auto& vss:vsss) r.push_back(to_vec3d2(vss)); return r;
}
static inline VectorPI1 to_pi1(const std::vector<std::pair<int,int>>& v) { return v; }
static inline std::vector<std::pair<int,int>> from_pi1(const VectorPI1& v) { return v; }
static inline VectorPB1 to_pb1(const std::vector<std::pair<bool,bool>>& v) { return v; }
// ============================================================
PYBIND11_MODULE(hgp_py, m) {
    m.doc() = "pybind11 bindings for libhgp — 完整版";

    //new
    py::class_<HGPMesh>(m, "HGPMesh")
    .def(py::init<>())
    .def_readwrite("verts",      &HGPMesh::verts)
    .def_readwrite("face_id_0",  &HGPMesh::face_id_0)
    .def_readwrite("face_id_1",  &HGPMesh::face_id_1)
    .def_readwrite("face_id_2",  &HGPMesh::face_id_2)
    .def("to_obj_string",        &HGPMesh::ToOBJString)
    .def("save_obj",             &HGPMesh::SaveOBJ)
    .def_static("load_obj",      &HGPMesh::LoadOBJ);
    //new end

    // ════════════════════════════════════════════════════
    // ── IO ───────────────────────────────────────────────
    // ════════════════════════════════════════════════════

    // HGP_Vector_Base: 求与 n 垂直的向量 r
    m.def("HGP_Vector_Base",
        [](const std::vector<double>& n) -> std::vector<double> {
            Vector3d r;
            HGP_Vector_Base(to_vec3d(n), r);
            return from_vec3d(r);
        }, py::arg("n"), "返回与输入向量垂直的单位向量 [x,y,z]。");

    // HGP_Test_PGL: 调试输出，原样暴露
    m.def("HGP_Test_PGL",
        [](const std::vector<double>& n, const std::string& s, const std::string& c){
            HGP_Test_PGL(to_vec3d(n), s.c_str(), c.c_str());
        }, py::arg("n"), py::arg("str"), py::arg("char_"), "调试函数，输出 PGL 信息。");

    // ════════════════════════════════════════════════════
    // ── HGP2D
    // ════════════════════════════════════════════════════
    m.def("HGP_2D_Square_Regular_Sampling_C1",
        [](double d) -> std::vector<std::vector<double>> {
            return from_vec2d1(HGP_2D_Square_Regular_Sampling_C1(d));
        }, py::arg("d"), "正方形域规则采样C1。");

    m.def("HGP_2D_Square_Regular_Sampling_C2",
        [](double d) -> std::pair<std::vector<std::vector<double>>,
                                   std::vector<std::pair<int,int>>> {
            VectorPI1 nb;
            auto pts = HGP_2D_Square_Regular_Sampling_C2(d, nb);
            return {from_vec2d1(pts), nb};
        }, py::arg("d"), "正方形域规则采样C2，返回 (采样点, 邻居对)。");

    m.def("HGP_2D_Square_Regular_Sampling_C3",
        [](double d, bool compute_neighbors)
           -> std::pair<std::vector<std::vector<double>>,
                        std::vector<std::pair<int,int>>> {
            VectorPI1 nb;
            auto pts = HGP_2D_Square_Regular_Sampling_C3(d, nb, compute_neighbors);
            return {from_vec2d1(pts), nb};
        }, py::arg("d"), py::arg("compute_neighbors"), "正方形域规则采样C3。");

    // HGP_2D_Distance_Point_Polygons
    m.def("HGP_2D_Distance_Point_Polygons",
        [](const std::vector<double>& p,
           const std::vector<std::vector<std::vector<double>>>& pys) -> double {
            return HGP_2D_Distance_Point_Polygons(to_vec2d(p), to_vec2d2(pys));
        }, py::arg("p"), py::arg("pys"), "2D点到多个多边形的最近距离。");

    // HGP_2D_Intersection_Segment_Line（线段与直线）
    m.def("HGP_2D_Intersection_Segment_Line",
        [](const std::vector<double>& ss, const std::vector<double>& se,
           const std::vector<double>& ls, const std::vector<double>& le)
           -> std::pair<bool, std::vector<double>> {
            Vector2d inter;
            bool hit = HGP_2D_Intersection_Segment_Line(
                to_vec2d(ss), to_vec2d(se), to_vec2d(ls), to_vec2d(le), inter);
            return {hit, from_vec2d(inter)};
        }, py::arg("s_s"), py::arg("s_e"), py::arg("l_s"), py::arg("l_e"),
        "2D线段与直线求交。返回 (是否相交, 交点[x,y])。");

    // HGP_2D_Intersection_Ray_Segment（射线与线段）
    m.def("HGP_2D_Intersection_Ray_Segment",
        [](const std::vector<double>& s0s, const std::vector<double>& s0e,
           const std::vector<double>& s1s, const std::vector<double>& s1e)
           -> std::pair<bool, std::vector<double>> {
            Vector2d inter;
            bool hit = HGP_2D_Intersection_Ray_Segment(
                to_vec2d(s0s), to_vec2d(s0e), to_vec2d(s1s), to_vec2d(s1e), inter);
            return {hit, from_vec2d(inter)};
        }, py::arg("s_0_s"), py::arg("s_0_e"), py::arg("s_1_s"), py::arg("s_1_e"),
        "2D射线与线段求交。返回 (是否相交, 交点[x,y])。");

    // HGP_2D_Intersection_Segment_Polygon
    m.def("HGP_2D_Intersection_Segment_Polygon",
        [](const std::vector<double>& ss, const std::vector<double>& se,
           const std::vector<std::vector<double>>& poly) -> bool {
            return HGP_2D_Intersection_Segment_Polygon(
                to_vec2d(ss), to_vec2d(se), to_vec2d1(poly));
        }, py::arg("s_s"), py::arg("s_e"), py::arg("p"), "2D线段与多边形是否相交。");

    // HGP_2D_Intersection_Polygon_Polygon
    m.def("HGP_2D_Intersection_Polygon_Polygon",
        [](const std::vector<std::vector<double>>& p1,
           const std::vector<std::vector<double>>& p2) -> bool {
            return HGP_2D_Intersection_Polygon_Polygon(to_vec2d1(p1), to_vec2d1(p2));
        }, py::arg("p1"), py::arg("p2"), "两个2D多边形是否相交。");

    // HGP_2D_Two_Polygons_Union（out: inter_polygons）
    m.def("HGP_2D_Two_Polygons_Union",
        [](const std::vector<std::vector<double>>& poly0,
           const std::vector<std::vector<double>>& poly1)
           -> std::pair<double, std::vector<std::vector<std::vector<double>>>> {
            Vector2d2 inter;
            double area = HGP_2D_Two_Polygons_Union(
                to_vec2d1(poly0), to_vec2d1(poly1), inter);
            return {area, from_vec2d2(inter)};
        }, py::arg("poly_0"), py::arg("poly_1"),
        "两多边形并集。返回 (并集面积, 并集多边形列表)。");

    // HGP_2D_Two_Polygons_Intersection（纯输入，返回面积）
    m.def("HGP_2D_Two_Polygons_Intersection",
        [](const std::vector<std::vector<double>>& poly0,
           const std::vector<std::vector<double>>& poly1) -> double {
            return HGP_2D_Two_Polygons_Intersection(to_vec2d1(poly0), to_vec2d1(poly1));
        }, py::arg("poly_0"), py::arg("poly_1"), "返回两多边形交集面积。");

    // HGP_Decompose_Polyline（out: result list[int]）
    m.def("HGP_Decompose_Polyline",
        [](const std::vector<std::vector<double>>& polyline, double threshold)
           -> std::vector<int> {
            Vector1i1 result;
            HGP_Decompose_Polyline(to_vec2d1(polyline), threshold, result);
            return result;
        }, py::arg("polyline"), py::arg("threshold"), "折线分解，返回关键点索引列表。");

    // HGP_Identify_Polycut_Extend（out: ns, ne）
    m.def("HGP_Identify_Polycut_Extend",
        [](const std::vector<std::vector<double>>& polygon,
           const std::vector<double>& s, const std::vector<double>& e)
           -> std::tuple<bool, std::vector<double>, std::vector<double>> {
            Vector2d ns, ne;
            bool ok = HGP_Identify_Polycut_Extend(
                to_vec2d1(polygon), to_vec2d(s), to_vec2d(e), ns, ne);
            return {ok, from_vec2d(ns), from_vec2d(ne)};
        }, py::arg("polygon"), py::arg("s"), py::arg("e"),
        "识别延伸的多边形切割线。返回 (bool, ns[x,y], ne[x,y])。");

    // HGP_Identify_Polycut_NotExtend（纯判断）
    m.def("HGP_Identify_Polycut_NotExtend",
        [](const std::vector<std::vector<double>>& polygon,
           const std::vector<double>& s, const std::vector<double>& e) -> bool {
            return HGP_Identify_Polycut_NotExtend(
                to_vec2d1(polygon), to_vec2d(s), to_vec2d(e));
        }, py::arg("polygon"), py::arg("s"), py::arg("e"),
        "判断切割线（非延伸）是否有效。");

    // HGP_Identify_Polycut（out: result VectorPB1）
    m.def("HGP_Identify_Polycut",
        [](const std::vector<std::vector<double>>& polygon,
           const std::vector<std::vector<double>>& cutLine)
           -> std::pair<bool, std::vector<std::pair<bool,bool>>> {
            VectorPB1 result;
            bool ok = HGP_Identify_Polycut(to_vec2d1(polygon), to_vec2d1(cutLine), result);
            return {ok, result};
        }, py::arg("polygon"), py::arg("cutLine"),
        "识别多边形切割。返回 (bool, list[(bool,bool)])。");

    // HGP_Construct_InOutSide_Polygon（out: isPInside, isQInside）
    m.def("HGP_Construct_InOutSide_Polygon",
        [](const std::vector<std::vector<double>>& py_,
           const std::vector<double>& p, const std::vector<double>& q)
           -> std::tuple<bool, bool, bool> {
            bool ip=false, iq=false;
            bool ok = HGP_Construct_InOutSide_Polygon(
                to_vec2d1(py_), to_vec2d(p), to_vec2d(q), ip, iq);
            return {ok, ip, iq};
        }, py::arg("py"), py::arg("p"), py::arg("q"),
        "判断 p、q 是否在多边形内。返回 (ok, isPInside, isQInside)。");

    // HGP_Get_Angle_Kerf_Offset_Tan
    m.def("HGP_Get_Angle_Kerf_Offset_Tan",
        [](const std::vector<double>& a, const std::vector<double>& b) -> double {
            return HGP_Get_Angle_Kerf_Offset_Tan(to_vec2d(a), to_vec2d(b));
        }, py::arg("a"), py::arg("b"), "计算两向量的 kerf offset 偏置角正切值。");

    // HGP_2D_Projection_Point_Segment
    m.def("HGP_2D_Projection_Point_Segment",
        [](const std::vector<double>& p,
           const std::vector<double>& s, const std::vector<double>& e)
           -> std::vector<double> {
            return from_vec2d(HGP_2D_Projection_Point_Segment(
                to_vec2d(p), to_vec2d(s), to_vec2d(e)));
        }, py::arg("p"), py::arg("s"), py::arg("e"),
        "2D点在线段上的投影点。");

    // HGP_2D_Detect_Polygon_Inside_C1/C2/C3/C4/C5
    m.def("HGP_2D_Detect_Polygon_Inside_C1",
        [](const std::vector<std::vector<double>>& outside_py,
           const std::vector<double>& p) -> bool {
            return HGP_2D_Detect_Polygon_Inside_C1(to_vec2d1(outside_py), to_vec2d(p));
        }, py::arg("outside_py"), py::arg("p"), "检测点 p 是否在多边形内（C1）。");

    m.def("HGP_2D_Detect_Polygon_Inside_C2",
        [](const std::vector<std::vector<double>>& outside_py,
           const std::vector<std::vector<double>>& inside_py) -> bool {
            return HGP_2D_Detect_Polygon_Inside_C2(
                to_vec2d1(outside_py), to_vec2d1(inside_py));
        }, py::arg("outside_py"), py::arg("inside_py"),
        "检测内多边形是否在外多边形内（C2）。");

    m.def("HGP_2D_Detect_Polygon_Inside_C3",
        [](const std::vector<std::vector<std::vector<double>>>& outside_pys,
           const std::vector<double>& p) -> bool {
            return HGP_2D_Detect_Polygon_Inside_C3(to_vec2d2(outside_pys), to_vec2d(p));
        }, py::arg("outside_pys"), py::arg("p"),
        "检测点 p 是否在多个外多边形内（C3）。");

    m.def("HGP_2D_Detect_Polygon_Inside_C4",
        [](const std::vector<std::vector<std::vector<double>>>& outside_pys,
           const std::vector<std::vector<double>>& inside_py) -> bool {
            return HGP_2D_Detect_Polygon_Inside_C4(
                to_vec2d2(outside_pys), to_vec2d1(inside_py));
        }, py::arg("outside_pys"), py::arg("inside_py"),
        "检测内多边形是否在多个外多边形内（C4）。");

    m.def("HGP_2D_Detect_Polygon_Inside_C5",
        [](const std::vector<std::vector<std::vector<double>>>& outside_pys,
           const std::vector<std::vector<std::vector<double>>>& inside_pys) -> bool {
            return HGP_2D_Detect_Polygon_Inside_C5(
                to_vec2d2(outside_pys), to_vec2d2(inside_pys));
        }, py::arg("outside_pys"), py::arg("inside_pys"),
        "检测多个内多边形是否在多个外多边形内（C5）。");

    // HGP_2D_Distance_Polygon_Polygon / Polygons_Polygons
    m.def("HGP_2D_Distance_Polygon_Polygon",
        [](const std::vector<std::vector<double>>& p0,
           const std::vector<std::vector<double>>& p1) -> double {
            return HGP_2D_Distance_Polygon_Polygon(to_vec2d1(p0), to_vec2d1(p1));
        }, py::arg("poly_0"), py::arg("poly_1"), "两个2D多边形之间的最近距离。");

    m.def("HGP_2D_Distance_Polygons_Polygons",
        [](const std::vector<std::vector<std::vector<double>>>& p0,
           const std::vector<std::vector<std::vector<double>>>& p1) -> double {
            return HGP_2D_Distance_Polygons_Polygons(to_vec2d2(p0), to_vec2d2(p1));
        }, py::arg("poly_0"), py::arg("poly_1"), "两组多边形之间的最近距离。");

    // HGP_2D_Nearest_Point_Polygon_C2（out: p, min_d）
    m.def("HGP_2D_Nearest_Point_Polygon_C2",
        [](const std::vector<double>& v,
           const std::vector<std::vector<double>>& poly)
           -> std::pair<std::vector<double>, double> {
            Vector2d p; double min_d=0.0;
            HGP_2D_Nearest_Point_Polygon_C2(to_vec2d(v), to_vec2d1(poly), p, min_d);
            return {from_vec2d(p), min_d};
        }, py::arg("v"), py::arg("poly"),
        "多边形上距 v 最近点（C2）。返回 (最近点[x,y], 最近距离)。");

    // HGP_2D_Nearest_Point_Polygons
    m.def("HGP_2D_Nearest_Point_Polygons",
        [](const std::vector<double>& v,
           const std::vector<std::vector<std::vector<double>>>& polys)
           -> std::vector<double> {
            return from_vec2d(HGP_2D_Nearest_Point_Polygons(
                to_vec2d(v), to_vec2d2(polys)));
        }, py::arg("v"), py::arg("polys"), "多组多边形中距 v 最近的点。");

    // HGP_2D_Polygon_Inside_Point_C1 / C2
    m.def("HGP_2D_Polygon_Inside_Point_C1",
        [](const std::vector<std::vector<double>>& poly) -> std::vector<double> {
            return from_vec2d(HGP_2D_Polygon_Inside_Point_C1(to_vec2d1(poly)));
        }, py::arg("poly"), "返回多边形内部的一个点（C1）。");

    m.def("HGP_2D_Polygon_Inside_Point_C2",
        [](const std::vector<std::vector<std::vector<double>>>& polys)
           -> std::pair<bool, std::vector<double>> {
            Vector2d inner;
            bool ok = HGP_2D_Polygon_Inside_Point_C2(to_vec2d2(polys), inner);
            return {ok, from_vec2d(inner)};
        }, py::arg("polys"),
        "返回多组多边形内部的一个点（C2）。返回 (bool, [x,y])。");

    // HGP_2D_Polygons_One_Offsets
    m.def("HGP_2D_Polygons_One_Offsets",
        [](const std::vector<std::vector<std::vector<double>>>& polys, double d)
           -> std::vector<std::vector<std::vector<double>>> {
            Vector2d2 out;
            HGP_2D_Polygons_One_Offsets(to_vec2d2(polys), d, out);
            return from_vec2d2(out);
        }, py::arg("polys"), py::arg("d"), "多个多边形的偏置（offset）。");

    // HGP_2D_Polygons_Simple / Polygon_Simple / Polygon_Simple_Inter
    m.def("HGP_2D_Polygons_Simple",
        [](const std::vector<std::vector<std::vector<double>>>& poly) -> bool {
            return HGP_2D_Polygons_Simple(to_vec2d2(poly));
        }, py::arg("poly"), "判断多组多边形是否均为简单多边形。");

    m.def("HGP_2D_Polygon_Simple",
        [](const std::vector<std::vector<double>>& poly) -> bool {
            return HGP_2D_Polygon_Simple(to_vec2d1(poly));
        }, py::arg("poly"), "判断单个多边形是否为简单多边形。");

    m.def("HGP_2D_Polygon_Simple_Inter",
        [](const std::vector<std::vector<double>>& poly) -> bool {
            return HGP_2D_Polygon_Simple_Inter(to_vec2d1(poly));
        }, py::arg("poly"), "判断多边形是否为简单多边形（基于自交检测）。");

    // HGP_2D_OBB_Box（out: center, axis_0, axis_1, extent_0, extent_1）
    m.def("HGP_2D_OBB_Box",
        [](const std::vector<std::vector<double>>& vec)
           -> std::tuple<std::vector<double>, std::vector<double>,
                         std::vector<double>, double, double> {
            Vector2d center, axis0, axis1; double e0=0, e1=0;
            HGP_2D_OBB_Box(to_vec2d1(vec), center, axis0, axis1, e0, e1);
            return {from_vec2d(center), from_vec2d(axis0), from_vec2d(axis1), e0, e1};
        }, py::arg("vec"),
        "计算2D点集的有向包围盒（OBB）。返回 (中心, 轴0, 轴1, 半长0, 半长1)。");

    // HGP_Image_Grid_Decomposition_C1（out: image, boundary_xs, boundary_ys）
    m.def("HGP_Image_Grid_Decomposition_C1",
        []() -> std::tuple<std::vector<std::vector<int>>,
                           std::vector<std::vector<double>>,
                           std::vector<std::vector<double>>> {
            Vector1i2 image; Vector1d2 bxs, bys;
            HGP_Image_Grid_Decomposition_C1(image, bxs, bys);
            return {image, bxs, bys};
        }, "图像网格分解C1。返回 (image, boundary_xs, boundary_ys)。");

    // HGP_Image_Grid_Decomposition_Conservative_C1
    m.def("HGP_Image_Grid_Decomposition_Conservative_C1",
        []() -> std::tuple<std::vector<std::vector<int>>,
                           std::vector<std::vector<double>>,
                           std::vector<std::vector<double>>> {
            Vector1i2 image; Vector1d2 bxs, bys;
            HGP_Image_Grid_Decomposition_Conservative_C1(image, bxs, bys);
            return {image, bxs, bys};
        }, "保守图像网格分解C1。返回 (image, boundary_xs, boundary_ys)。");

    // HGP_Image_Grid_Decomposition_C2
    m.def("HGP_Image_Grid_Decomposition_C2",
        []() -> std::pair<std::vector<std::vector<int>>,
                          std::vector<std::vector<std::vector<double>>>> {
            Vector1i2 image; Vector2d2 boundaries;
            HGP_Image_Grid_Decomposition_C2(image, boundaries);
            return {image, from_vec2d2(boundaries)};
        }, "图像网格分解C2。返回 (image, boundaries)。");

    // HGP_Image_Grid_Decomposition_Conservative_C2
    m.def("HGP_Image_Grid_Decomposition_Conservative_C2",
        []() -> std::pair<std::vector<std::vector<int>>,
                          std::vector<std::vector<std::vector<double>>>> {
            Vector1i2 image; Vector2d2 boundaries;
            HGP_Image_Grid_Decomposition_Conservative_C2(image, boundaries);
            return {image, from_vec2d2(boundaries)};
        }, "保守图像网格分解C2。返回 (image, boundaries)。");

    // ════════════════════════════════════════════════════
    // ── HGP3D
    // ════════════════════════════════════════════════════

    // HGP_3D_Plane_Point_Projection（out: result）
    m.def("HGP_3D_Plane_Point_Projection",
        [](const std::vector<double>& pp, const std::vector<double>& pn,
           const std::vector<double>& p) -> std::vector<double> {
            Vector3d result;
            HGP_3D_Plane_Point_Projection(to_vec3d(pp), to_vec3d(pn), to_vec3d(p), result);
            return from_vec3d(result);
        }, py::arg("plane_p"), py::arg("plane_n"), py::arg("p"),
        "3D点投影到平面，返回投影点 [x,y,z]。");

    // HGP_3D_Plane_Points_Projection（out: project_points）
    m.def("HGP_3D_Plane_Points_Projection",
        [](const std::vector<double>& pp, const std::vector<double>& pn,
           const std::vector<std::vector<double>>& pts)
           -> std::vector<std::vector<double>> {
            Vector3d1 out;
            HGP_3D_Plane_Points_Projection(to_vec3d(pp), to_vec3d(pn), to_vec3d1(pts), out);
            return from_vec3d1(out);
        }, py::arg("plane_p"), py::arg("plane_n"), py::arg("points"),
        "3D点集批量投影到平面。");

    // HGP_3D_Plane_3D_to_2D_Point（out: result Vector2d）
    m.def("HGP_3D_Plane_3D_to_2D_Point",
        [](const std::vector<double>& pp, const std::vector<double>& pn,
           const std::vector<double>& pt3d) -> std::vector<double> {
            Vector2d result;
            HGP_3D_Plane_3D_to_2D_Point(to_vec3d(pp), to_vec3d(pn), to_vec3d(pt3d), result);
            return from_vec2d(result);
        }, py::arg("plane_p"), py::arg("plane_n"), py::arg("point_3d"),
        "3D点投影到平面局部2D坐标。");

    // HGP_3D_Plane_2D_to_3D_Point（out: result Vector3d）
    m.def("HGP_3D_Plane_2D_to_3D_Point",
        [](const std::vector<double>& pp, const std::vector<double>& pn,
           const std::vector<double>& pt2d) -> std::vector<double> {
            Vector3d result;
            HGP_3D_Plane_2D_to_3D_Point(to_vec3d(pp), to_vec3d(pn), to_vec2d(pt2d), result);
            return from_vec3d(result);
        }, py::arg("plane_p"), py::arg("plane_n"), py::arg("points_2d"),
        "平面局部2D坐标反投影��3D点。");

    // HGP_3D_Plane_3D_to_2D_Points（batch, out: points_2d）
    m.def("HGP_3D_Plane_3D_to_2D_Points",
        [](const std::vector<double>& pp, const std::vector<double>& pn,
           const std::vector<std::vector<double>>& pts3d)
           -> std::vector<std::vector<double>> {
            Vector2d1 out;
            HGP_3D_Plane_3D_to_2D_Points(to_vec3d(pp), to_vec3d(pn), to_vec3d1(pts3d), out);
            return from_vec2d1(out);
        }, py::arg("plane_p"), py::arg("plane_n"), py::arg("points_3d"),
        "批量3D点投影为平面2D坐标。");

    // HGP_3D_Plane_3Ds_to_2Ds_Points（batch batch, out: points_2d）
    m.def("HGP_3D_Plane_3Ds_to_2Ds_Points",
        [](const std::vector<double>& pp, const std::vector<double>& pn,
           const std::vector<std::vector<std::vector<double>>>& pts3d)
           -> std::vector<std::vector<std::vector<double>>> {
            Vector2d2 out;
            HGP_3D_Plane_3Ds_to_2Ds_Points(to_vec3d(pp), to_vec3d(pn), to_vec3d2(pts3d), out);
            return from_vec2d2(out);
        }, py::arg("plane_p"), py::arg("plane_n"), py::arg("points_3d"),
        "批量多组3D点���影为平面2D坐标。");

    // HGP_3D_Plane_2D_to_3D_Points（batch, out: points_3d）
    m.def("HGP_3D_Plane_2D_to_3D_Points",
        [](const std::vector<double>& pp, const std::vector<double>& pn,
           const std::vector<std::vector<double>>& pts2d)
           -> std::vector<std::vector<double>> {
            Vector3d1 out;
            HGP_3D_Plane_2D_to_3D_Points(to_vec3d(pp), to_vec3d(pn), to_vec2d1(pts2d), out);
            return from_vec3d1(out);
        }, py::arg("plane_p"), py::arg("plane_n"), py::arg("points_2d"),
        "批量平面2D坐标反投影为3D点。");

    // HGP_3D_Projection_Point_Segment
    m.def("HGP_3D_Projection_Point_Segment",
        [](const std::vector<double>& p,
           const std::vector<double>& ss, const std::vector<double>& se)
           -> std::vector<double> {
            return from_vec3d(HGP_3D_Projection_Point_Segment(
                to_vec3d(p), to_vec3d(ss), to_vec3d(se)));
        }, py::arg("p"), py::arg("s_s"), py::arg("s_e"),
        "3D点在线段上的投影点。");

    // HGP_3D_Distance_Point_Polygon（注意参数顺序是 py, p）
    m.def("HGP_3D_Distance_Point_Polygon",
        [](const std::vector<std::vector<double>>& py_,
           const std::vector<double>& p) -> double {
            return HGP_3D_Distance_Point_Polygon(to_vec3d1(py_), to_vec3d(p));
        }, py::arg("py"), py::arg("p"), "3D点到多边形的距离（注意参数顺序：先多边形后点）。");

    // HGP_2D_Polygon_Triangulation_C1（out: faces Vector1i2）
    m.def("HGP_2D_Polygon_Triangulation_C1",
        [](const std::vector<std::vector<std::vector<double>>>& polys)
           -> std::vector<std::vector<int>> {
            Vector1i2 faces;
            HGP_2D_Polygon_Triangulation_C1(to_vec2d2(polys), faces);
            return faces;
        }, py::arg("polys"),
        "2D多边形三角剖分C1。返回三角面索引列表 list[list[int]]（每行3个索引）。");

    // HGP_2D_Polygon_Triangulation_C2（直接返回 Vector1i2）
    m.def("HGP_2D_Polygon_Triangulation_C2",
        [](const std::vector<std::vector<std::vector<double>>>& polys)
           -> std::vector<std::vector<int>> {
            return HGP_2D_Polygon_Triangulation_C2(to_vec2d2(polys));
        }, py::arg("polys"), "2D多边形三角剖分C2（直接返回）。");

    // HGP_2D_Polygon_Triangulation_C3（单多边形）
    m.def("HGP_2D_Polygon_Triangulation_C3",
        [](const std::vector<std::vector<double>>& poly)
           -> std::vector<std::vector<int>> {
            return HGP_2D_Polygon_Triangulation_C3(to_vec2d1(poly));
        }, py::arg("poly"), "单个2D多边形三角剖分C3。");

    // HGP_3D_Distance_Point_Line
    m.def("HGP_3D_Distance_Point_Line",
        [](const std::vector<double>& p,
           const std::vector<double>& ls, const std::vector<double>& le) -> double {
            return HGP_3D_Distance_Point_Line(to_vec3d(p), to_vec3d(ls), to_vec3d(le));
        }, py::arg("p"), py::arg("l_s"), py::arg("l_e"), "3D点到直线的距离。");

    // HGP_3D_Projection_Point_Line
    m.def("HGP_3D_Projection_Point_Line",
        [](const std::vector<double>& p,
           const std::vector<double>& ls, const std::vector<double>& le)
           -> std::vector<double> {
            return from_vec3d(HGP_3D_Projection_Point_Line(
                to_vec3d(p), to_vec3d(ls), to_vec3d(le)));
        }, py::arg("p"), py::arg("l_s"), py::arg("l_e"), "3D点在直线上的投影点。");

    // HGP_3D_Distance_Segment_Segment
    m.def("HGP_3D_Distance_Segment_Segment",
        [](const std::vector<double>& s0s, const std::vector<double>& s0e,
           const std::vector<double>& s1s, const std::vector<double>& s1e) -> double {
            return HGP_3D_Distance_Segment_Segment(
                to_vec3d(s0s), to_vec3d(s0e), to_vec3d(s1s), to_vec3d(s1e));
        }, py::arg("s_0_s"), py::arg("s_0_e"), py::arg("s_1_s"), py::arg("s_1_e"),
        "3D两线段之间的最近距离。");

    // HGP_3D_Distance_Point_Plane
    m.def("HGP_3D_Distance_Point_Plane",
        [](const std::vector<double>& v,
           const std::vector<double>& pp, const std::vector<double>& pn) -> double {
            return HGP_3D_Distance_Point_Plane(to_vec3d(v), to_vec3d(pp), to_vec3d(pn));
        }, py::arg("v"), py::arg("plane_p"), py::arg("plane_n"), "3D点到平面的距离。");

    // HGP_3D_Intersection_Segment_Line（out: inter）
    m.def("HGP_3D_Intersection_Segment_Line",
        [](const std::vector<double>& ss, const std::vector<double>& se,
           const std::vector<double>& ls, const std::vector<double>& le)
           -> std::pair<bool, std::vector<double>> {
            Vector3d inter;
            bool hit = HGP_3D_Intersection_Segment_Line(
                to_vec3d(ss), to_vec3d(se), to_vec3d(ls), to_vec3d(le), inter);
            return {hit, from_vec3d(inter)};
        }, py::arg("s_s"), py::arg("s_e"), py::arg("l_s"), py::arg("l_e"),
        "3D线段与直线求交。返回 (bool, [x,y,z])。");

    // HGP_3D_Intersection_Segment_Segment（out: iter）
    m.def("HGP_3D_Intersection_Segment_Segment",
        [](const std::vector<double>& s0s, const std::vector<double>& s0e,
           const std::vector<double>& s1s, const std::vector<double>& s1e)
           -> std::pair<bool, std::vector<double>> {
            Vector3d inter;
            bool hit = HGP_3D_Intersection_Segment_Segment(
                to_vec3d(s0s), to_vec3d(s0e), to_vec3d(s1s), to_vec3d(s1e), inter);
            return {hit, from_vec3d(inter)};
        }, py::arg("s_0_s"), py::arg("s_0_e"), py::arg("s_1_s"), py::arg("s_1_e"),
        "3D两线段求交。返回 (bool, [x,y,z])。");

    // HGP_3D_Intersection_Line_Plane（out: inter）
    m.def("HGP_3D_Intersection_Line_Plane",
        [](const std::vector<double>& ls, const std::vector<double>& le,
           const std::vector<double>& pp, const std::vector<double>& pn)
           -> std::pair<bool, std::vector<double>> {
            Vector3d inter;
            bool hit = HGP_3D_Intersection_Line_Plane(
                to_vec3d(ls), to_vec3d(le), to_vec3d(pp), to_vec3d(pn), inter);
            return {hit, from_vec3d(inter)};
        }, py::arg("l_s"), py::arg("l_e"), py::arg("plane_p"), py::arg("plane_n"),
        "3D直线与平面求交。返回 (bool, [x,y,z])。");

    // HGP_3D_Projection_Point_Plane_C2（三点定平面）
    m.def("HGP_3D_Projection_Point_Plane_C2",
        [](const std::vector<double>& p,
           const std::vector<double>& p0,
           const std::vector<double>& p1,
           const std::vector<double>& p2) -> std::vector<double> {
            return from_vec3d(HGP_3D_Projection_Point_Plane_C2(
                to_vec3d(p), to_vec3d(p0), to_vec3d(p1), to_vec3d(p2)));
        }, py::arg("p"), py::arg("plane_p_0"), py::arg("plane_p_1"), py::arg("plane_p_2"),
        "3D点投影到三点定义的平面（C2）。");

    // HGP_3D_Projection_3D_Point_Plane_2D_C1 / C2
    m.def("HGP_3D_Projection_3D_Point_Plane_2D_C1",
        [](const std::vector<double>& p,
           const std::vector<double>& pp, const std::vector<double>& pn)
           -> std::vector<double> {
            return from_vec2d(HGP_3D_Projection_3D_Point_Plane_2D_C1(
                to_vec3d(p), to_vec3d(pp), to_vec3d(pn)));
        }, py::arg("p"), py::arg("plane_p"), py::arg("plane_n"),
        "3D点投影到平面并输出局部2D坐标（C1）。");

    m.def("HGP_3D_Projection_3D_Point_Plane_2D_C2",
        [](const std::vector<double>& p,
           const std::vector<double>& p0,
           const std::vector<double>& p1,
           const std::vector<double>& p2) -> std::vector<double> {
            return from_vec2d(HGP_3D_Projection_3D_Point_Plane_2D_C2(
                to_vec3d(p), to_vec3d(p0), to_vec3d(p1), to_vec3d(p2)));
        }, py::arg("p"), py::arg("plane_p_0"), py::arg("plane_p_1"), py::arg("plane_p_2"),
        "3D点投影到平面并输出局部2D坐标（C2，三点定平面）。");

    // HGP_3D_Plane_ABCD（out: a,b,c,d）
    m.def("HGP_3D_Plane_ABCD",
        [](const std::vector<double>& pp, const std::vector<double>& pn)
           -> std::tuple<double,double,double,double> {
            double a=0,b=0,c=0,d=0;
            HGP_3D_Plane_ABCD(to_vec3d(pp), to_vec3d(pn), a, b, c, d);
            return {a,b,c,d};
        }, py::arg("plane_p"), py::arg("plane_n"),
        "计算平面方程 ax+by+cz+d=0 的系数。返回 (a,b,c,d)。");

    // HGP_3D_Plane_Base_1
    m.def("HGP_3D_Plane_Base_1",
        [](const std::vector<double>& pp, const std::vector<double>& pn)
           -> std::vector<double> {
            return from_vec3d(HGP_3D_Plane_Base_1(to_vec3d(pp), to_vec3d(pn)));
        }, py::arg("plane_p"), py::arg("plane_n"),
        "返回平面的一个正交基向量（在平面内，垂直于法向量）。");

    // HGP_Face_Normal（out: normal_0/1/2，返回值 + 3个out）
    m.def("HGP_Face_Normal",
        [](const std::vector<double>& source,
           const std::vector<double>& t0,
           const std::vector<double>& t1,
           const std::vector<double>& t2)
           -> std::tuple<std::vector<double>,
                         std::vector<double>,
                         std::vector<double>,
                         std::vector<double>> {
            Vector3d n0, n1, n2;
            Vector3d ret = HGP_Face_Normal(
                to_vec3d(source), to_vec3d(t0), to_vec3d(t1), to_vec3d(t2), n0, n1, n2);
            return {from_vec3d(ret), from_vec3d(n0), from_vec3d(n1), from_vec3d(n2)};
        }, py::arg("source"), py::arg("tri_0"), py::arg("tri_1"), py::arg("tri_2"),
        "计算三角面法向量。返回 (face_normal, normal_0, normal_1, normal_2)。");

    // ════════════════════════════════════════════════════
    // ── HGPMESH
    // ════════════════════════════════════════════════════

    // HGP_Remesh_Surface_by_Adding_Feature
    // out: igl_cutting_0_edges, igl_cutting_1_edges, igl_cutting_points, cutting_faces
    m.def("HGP_Remesh_Surface_by_Adding_Feature",
        [](const std::vector<std::vector<double>>& feature,
           const std::vector<int>& face_ids,
           const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0,
           const std::vector<int>& f1,
           const std::vector<int>& f2)
           -> std::tuple<std::vector<int>, std::vector<int>,
                         std::vector<std::vector<double>>,
                         std::vector<std::vector<int>>> {
            Vector1i1 e0, e1; Vector3d1 pts; Vector1i2 cfaces;
            HGP_Remesh_Surface_by_Adding_Feature(
                to_vec3d1(feature), face_ids, to_vec3d1(vecs),
                f0, f1, f2, e0, e1, pts, cfaces);
            return {e0, e1, from_vec3d1(pts), cfaces};
        }, py::arg("feature"), py::arg("face_ids"), py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "网格特征线重网格化。返回 (edge0, edge1, cutting_points, cutting_faces)。");

    // HGP_Mesh_Edges（仅 path，无返回值，副作用写文件）
    m.def("HGP_Mesh_Edges",
        [](const std::string& path){ HGP_Mesh_Edges(path.c_str()); },
        py::arg("path"), "计算并输出网格的边信息（写文件）。");

    // HGP_3D_Intersection_Sphere_Ray（out: i_x, i_y, i_z）
    m.def("HGP_3D_Intersection_Sphere_Ray",
        [](double cx, double cy, double cz, double radius,
           double ox, double oy, double oz,
           double dx, double dy, double dz)
           -> std::tuple<bool, std::vector<double>,
                         std::vector<double>, std::vector<double>> {
            Vector1d1 ix, iy, iz;
            bool hit = HGP_3D_Intersection_Sphere_Ray(
                cx,cy,cz,radius, ox,oy,oz, dx,dy,dz, ix,iy,iz);
            return {hit, ix, iy, iz};
        }, py::arg("center_x"), py::arg("center_y"), py::arg("center_z"),
           py::arg("radius"),
           py::arg("ray_origin_x"), py::arg("ray_origin_y"), py::arg("ray_origin_z"),
           py::arg("ray_direction_x"), py::arg("ray_direction_y"), py::arg("ray_direction_z"),
        "球与射线求交。返回 (bool, i_x[], i_y[], i_z[]) 各交点坐标分量列表。");

    // HGP_3D_Intersection_Ray_Triangle
    m.def("HGP_3D_Intersection_Ray_Triangle",
        [](const std::vector<double>& p, const std::vector<double>& n,
           const std::vector<double>& p0,
           const std::vector<double>& p1,
           const std::vector<double>& p2) -> bool {
            return HGP_3D_Intersection_Ray_Triangle(
                to_vec3d(p), to_vec3d(n),
                to_vec3d(p0), to_vec3d(p1), to_vec3d(p2));
        }, py::arg("p"), py::arg("n"),
           py::arg("p0"), py::arg("p1"), py::arg("p2"),
        "射线与三角形求交（布尔）。");

    // HGP_3D_Intersection_Ray_Mesh / Segment_Mesh
    m.def("HGP_3D_Intersection_Ray_Mesh",
        [](const std::vector<double>& p, const std::vector<double>& n,
           const std::string& path) -> bool {
            return HGP_3D_Intersection_Ray_Mesh(to_vec3d(p), to_vec3d(n), path.c_str());
        }, py::arg("p"), py::arg("n"), py::arg("path"), "射线与网格文件求交（布尔）。");

    m.def("HGP_3D_Intersection_Segment_Mesh",
        [](const std::vector<double>& s, const std::vector<double>& e,
           const std::string& path) -> bool {
            return HGP_3D_Intersection_Segment_Mesh(to_vec3d(s), to_vec3d(e), path.c_str());
        }, py::arg("s"), py::arg("e"), py::arg("path"), "线段与网格文件求交（布尔）。");

    // HGP_3D_Intersection_Segments_Mesh（out: inters Vector1b1）
    m.def("HGP_3D_Intersection_Segments_Mesh",
        [](const std::vector<std::vector<double>>& ss,
           const std::vector<std::vector<double>>& ee,
           const std::string& path) -> std::vector<bool> {
            Vector1b1 inters;
            HGP_3D_Intersection_Segments_Mesh(to_vec3d1(ss), to_vec3d1(ee), path.c_str(), inters);
            return inters;
        }, py::arg("ss"), py::arg("ee"), py::arg("path"),
        "批量线段与网格求交。返回 bool 列表。");

    // HGP_3D_Intersection_Polygons_Mesh（out: inters）
    m.def("HGP_3D_Intersection_Polygons_Mesh",
        [](const std::vector<std::vector<std::vector<double>>>& polygons,
           const std::string& path) -> std::vector<bool> {
            Vector1b1 inters;
            HGP_3D_Intersection_Polygons_Mesh(to_vec3d2(polygons), path.c_str(), inters);
            return inters;
        }, py::arg("polygons"), py::arg("path"),
        "批量多边形与网格求交。返回每个多边形是否相交的 bool 列表。");

    // HGP_3D_Intersection_Polygons_Mesh_Bool
    m.def("HGP_3D_Intersection_Polygons_Mesh_Bool",
        [](const std::vector<std::vector<std::vector<double>>>& polygons,
           const std::string& path) -> bool {
            return HGP_3D_Intersection_Polygons_Mesh_Bool(to_vec3d2(polygons), path.c_str());
        }, py::arg("polygons"), py::arg("path"),
        "是否有任意一个多边形与网格相交（布尔）。");

    // HGP_3D_Intersection_Rays_Mesh_Vector3d（out: inters Vector3d1）
    m.def("HGP_3D_Intersection_Rays_Mesh_Vector3d",
        [](const std::vector<std::vector<double>>& ps,
           const std::vector<std::vector<double>>& ns,
           const std::string& path) -> std::vector<std::vector<double>> {
            Vector3d1 inters;
            HGP_3D_Intersection_Rays_Mesh_Vector3d(
                to_vec3d1(ps), to_vec3d1(ns), path.c_str(), inters);
            return from_vec3d1(inters);
        }, py::arg("ps"), py::arg("ns"), py::arg("path"),
        "批量射线与网格求交，返回交点列表。");

    // HGP_3D_Intersection_Rays_Mesh_C1_Bool（每个点各自方向组, out: inters Vector1b2）
    m.def("HGP_3D_Intersection_Rays_Mesh_C1_Bool",
        [](const std::vector<std::vector<double>>& ps,
           const std::vector<std::vector<std::vector<double>>>& nes,
           const std::string& path) -> std::vector<std::vector<bool>> {
            Vector1b2 inters;
            HGP_3D_Intersection_Rays_Mesh_C1_Bool(
                to_vec3d1(ps), to_vec3d2(nes), path.c_str(), inters);
            return inters;
        }, py::arg("ps"), py::arg("nes"), py::arg("path"),
        "每个点用各自方向组测试（C1）。返回 list[list[bool]]。");

    // HGP_3D_Intersection_Rays_Mesh_C2_Bool（所有点共用方向集, out: inters Vector1b2）
    m.def("HGP_3D_Intersection_Rays_Mesh_C2_Bool",
        [](const std::vector<std::vector<double>>& ps,
           const std::vector<std::vector<double>>& ns,
           const std::string& path) -> std::vector<std::vector<bool>> {
            Vector1b2 inters;
            HGP_3D_Intersection_Rays_Mesh_C2_Bool(
                to_vec3d1(ps), to_vec3d1(ns), path.c_str(), inters);
            return inters;
        }, py::arg("ps"), py::arg("ns"), py::arg("path"),
        "所有点共用方向集测试（C2）。返回 list[list[bool]]。");

    // HGP_3D_Intersection_Rays_Mesh_C2_Vector3d（out: inters Vector1d2）
    m.def("HGP_3D_Intersection_Rays_Mesh_C2_Vector3d",
        [](const std::vector<std::vector<double>>& ps,
           const std::vector<std::vector<double>>& ns,
           const std::string& path) -> std::vector<std::vector<double>> {
            Vector1d2 inters;
            HGP_3D_Intersection_Rays_Mesh_C2_Vector3d(
                to_vec3d1(ps), to_vec3d1(ns), path.c_str(), inters);
            return inters;
        }, py::arg("ps"), py::arg("ns"), py::arg("path"),
        "所有点共用方向集，返回交点距离 list[list[float]]。");

    // HGP_3D_Points_Inside_Triangles_C1_Bool / C2_Bool
    m.def("HGP_3D_Points_Inside_Triangles_C1_Bool",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2,
           const std::vector<std::vector<double>>& points) -> std::vector<bool> {
            Vector1b1 insides;
            HGP_3D_Points_Inside_Triangles_C1_Bool(
                to_vec3d1(vecs), f0, f1, f2, to_vec3d1(points), insides);
            return insides;
        }, py::arg("vecs"), py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
           py::arg("points"),
        "判断点集是否在三角网格内（C1，数组传参）。");

    m.def("HGP_3D_Points_Inside_Triangles_C2_Bool",
        [](const std::string& path,
           const std::vector<std::vector<double>>& points) -> std::vector<bool> {
            Vector1b1 insides;
            HGP_3D_Points_Inside_Triangles_C2_Bool(path.c_str(), to_vec3d1(points), insides);
            return insides;
        }, py::arg("path"), py::arg("points"),
        "判断点集是否在网格内（C2，文件路径传参）。");

    // HGP_3D_Mesh_Dart_Sampling_C1 / C2（out: sampling_points）
    m.def("HGP_3D_Mesh_Dart_Sampling_C1",
        [](const std::string& outside_path, double d, int total_iter)
           -> std::vector<std::vector<double>> {
            Vector3d1 pts;
            HGP_3D_Mesh_Dart_Sampling_C1(outside_path.c_str(), d, pts, total_iter);
            return from_vec3d1(pts);
        }, py::arg("outside_path"), py::arg("d"), py::arg("total_iter"),
        "3D网格 dart 采样（C1，单层）。");

    m.def("HGP_3D_Mesh_Dart_Sampling_C2",
        [](const std::string& outside_path, const std::string& inside_path,
           double d, int total_iter) -> std::vector<std::vector<double>> {
            Vector3d1 pts;
            HGP_3D_Mesh_Dart_Sampling_C2(
                outside_path.c_str(), inside_path.c_str(), d, pts, total_iter);
            return from_vec3d1(pts);
        }, py::arg("outside_path"), py::arg("inside_path"),
           py::arg("d"), py::arg("total_iter"),
        "3D网格 dart 采样（C2，内外层）。");

    // HGP_3D_Mesh_Regular_Sampling_C1 / C2 / C3
    m.def("HGP_3D_Mesh_Regular_Sampling_C1",
        [](const std::string& outside_path, double d) -> std::vector<std::vector<double>> {
            Vector3d1 pts;
            HGP_3D_Mesh_Regular_Sampling_C1(outside_path.c_str(), d, pts);
            return from_vec3d1(pts);
        }, py::arg("outside_path"), py::arg("d"), "3D网格规则采样C1。");

    m.def("HGP_3D_Mesh_Regular_Sampling_C2",
        [](const std::string& outside_path, const std::string& inside_path, double d)
           -> std::vector<std::vector<double>> {
            Vector3d1 pts;
            HGP_3D_Mesh_Regular_Sampling_C2(outside_path.c_str(), inside_path.c_str(), d, pts);
            return from_vec3d1(pts);
        }, py::arg("outside_path"), py::arg("inside_path"), py::arg("d"),
        "3D网格规则采样C2（内外层）。");

    m.def("HGP_3D_Mesh_Regular_Sampling_C3",
        [](const std::string& outside_path, const std::string& inside_path, double d)
           -> std::pair<std::vector<std::vector<double>>,
                        std::vector<std::pair<int,int>>> {
            Vector3d1 pts; VectorPI1 nb;
            HGP_3D_Mesh_Regular_Sampling_C3(
                outside_path.c_str(), inside_path.c_str(), d, pts, nb);
            return {from_vec3d1(pts), nb};
        }, py::arg("outside_path"), py::arg("inside_path"), py::arg("d"),
        "3D网格规则采样C3（含邻居）。返回 (采样点, 邻居对)。");

    // HGP_3D_Cube_Surface_Sampling_C1 / C2 / C3
    // out: sampling_points (Vector3d2), neighbors (VectorPI2)
    m.def("HGP_3D_Cube_Surface_Sampling_C1",
        [](double cube_size, double d, bool compute_neighbors)
           -> std::pair<std::vector<std::vector<std::vector<double>>>,
                        std::vector<std::vector<std::pair<int,int>>>> {
            Vector3d2 pts; VectorPI2 nb;
            HGP_3D_Cube_Surface_Sampling_C1(cube_size, d, pts, nb, compute_neighbors);
            return {from_vec3d2(pts), nb};
        }, py::arg("cube_size"), py::arg("d"), py::arg("compute_neighbors"),
        "立方体表面采样C1。返回 (每面采样点列表, 邻居对列表)。");

    m.def("HGP_3D_Cube_Surface_Sampling_C2",
        [](double cube_size, double d)
           -> std::vector<std::vector<std::vector<double>>> {
            Vector3d2 pts;
            HGP_3D_Cube_Surface_Sampling_C2(cube_size, d, pts);
            return from_vec3d2(pts);
        }, py::arg("cube_size"), py::arg("d"), "立方体表面采样C2（不含邻居）。");

    m.def("HGP_3D_Cube_Surface_Sampling_C3",
        [](double cube_size, double d)
           -> std::pair<std::vector<std::vector<std::vector<double>>>,
                        std::vector<std::vector<std::pair<int,int>>>> {
            Vector3d2 pts; VectorPI2 nb;
            HGP_3D_Cube_Surface_Sampling_C3(cube_size, d, pts, nb);
            return {from_vec3d2(pts), nb};
        }, py::arg("cube_size"), py::arg("d"), "立方体表面采样C3（含邻居）。");

    // HGP_3D_Distance_Point_Triangle / Triangles
    m.def("HGP_3D_Distance_Point_Triangle",
        [](const std::vector<double>& p,
           const std::vector<double>& t0,
           const std::vector<double>& t1,
           const std::vector<double>& t2) -> double {
            return HGP_3D_Distance_Point_Triangle(
                to_vec3d(p), to_vec3d(t0), to_vec3d(t1), to_vec3d(t2));
        }, py::arg("p"), py::arg("t_0"), py::arg("t_1"), py::arg("t_2"),
        "3D点到三角形的距离。");

    m.def("HGP_3D_Distance_Point_Triangles",
        [](const std::vector<double>& p,
           const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0,
           const std::vector<int>& f1,
           const std::vector<int>& f2) -> double {
            return HGP_3D_Distance_Point_Triangles(
                to_vec3d(p), to_vec3d1(vecs), f0, f1, f2);
        }, py::arg("p"), py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "3D点到三角网格的最近距离。");

    // HGP_3D_Nearest_Point_Triangles
    m.def("HGP_3D_Nearest_Point_Triangles",
        [](const std::vector<double>& p,
           const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0,
           const std::vector<int>& f1,
           const std::vector<int>& f2) -> std::vector<double> {
            return from_vec3d(HGP_3D_Nearest_Point_Triangles(
                to_vec3d(p), to_vec3d1(vecs), f0, f1, f2));
        }, py::arg("p"), py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "3D点在三角网格上的最近点。");

    // HGP_3D_Distance_Point_Mesh / Neareast_Point_Mesh
    m.def("HGP_3D_Distance_Point_Mesh",
        [](const std::string& path,
           const std::vector<std::vector<double>>& qpts) -> std::vector<double> {
            Vector1d1 dists;
            HGP_3D_Distance_Point_Mesh(path.c_str(), to_vec3d1(qpts), dists);
            return dists;
        }, py::arg("path"), py::arg("query_points"),
        "查询点到网格的最近距离列表。");

    m.def("HGP_3D_Neareast_Point_Mesh",
        [](const std::string& path,
           const std::vector<std::vector<double>>& ves)
           -> std::vector<std::vector<double>> {
            Vector3d1 ners;
            HGP_3D_Neareast_Point_Mesh(path.c_str(), to_vec3d1(ves), ners);
            return from_vec3d1(ners);
        }, py::arg("path"), py::arg("ves"), "查询点在网格上的最近点列表。");

    // HGP_3D_Mesh_Near_Triangles（out: triangles Vector1i2）
    m.def("HGP_3D_Mesh_Near_Triangles",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2,
           const std::vector<std::vector<double>>& points, double d)
           -> std::vector<std::vector<int>> {
            Vector1i2 triangles;
            HGP_3D_Mesh_Near_Triangles(
                to_vec3d1(vecs), f0, f1, f2, to_vec3d1(points), d, triangles);
            return triangles;
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
           py::arg("points"), py::arg("d"),
        "查找每个点附近距离 d 内的三角形索引列表。");

    // HGP_3D_Points_inside_Triangles_C1 / C2（out: insides）
    m.def("HGP_3D_Points_inside_Triangles_C1",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2,
           const std::vector<std::vector<double>>& points) -> std::vector<bool> {
            Vector1b1 insides;
            HGP_3D_Points_inside_Triangles_C1(
                to_vec3d1(vecs), f0, f1, f2, to_vec3d1(points), insides);
            return insides;
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
           py::arg("points"), "判断点集是否在三角网格内（C1）。");

    m.def("HGP_3D_Points_inside_Triangles_C2",
        [](const std::string& path,
           const std::vector<std::vector<double>>& points) -> std::vector<bool> {
            Vector1b1 insides;
            HGP_3D_Points_inside_Triangles_C2(path.c_str(), to_vec3d1(points), insides);
            return insides;
        }, py::arg("path"), py::arg("points"), "判断点集是否在网格内（C2，文件）。");

    // HGP_Mesh_Subdivision / Loop_Subdivision_One_Step
    m.def("HGP_Mesh_Subdivision",
        [](const std::string& in, const std::string& sub, int step, const std::string& out){
            HGP_Mesh_Subdivision(in.c_str(), sub.c_str(), step, out.c_str());
        }, py::arg("in_path"), py::arg("sub"), py::arg("step"), py::arg("out_path"),
        "网格细分（文件IO）。sub: 细分类型字符串。");

    m.def("HGP_Mesh_Loop_Subdivision_One_Step",
        [](std::vector<std::vector<double>> vecs,
           std::vector<int> f0, std::vector<int> f1, std::vector<int> f2)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>> {
            Vector3d1 v = to_vec3d1(vecs);
            HGP_Mesh_Loop_Subdivision_One_Step(v, f0, f1, f2);
            return {from_vec3d1(v), f0, f1, f2};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "原地 Loop 细分一步。返回 (新顶点, 新face_id_0, 新face_id_1, 新face_id_2)。");

    // HGP_3D_Mesh_Curvature_C2（face_ids 用 Vector1i2）
    m.def("HGP_3D_Mesh_Curvature_C2",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<std::vector<int>>& face_ids)
           -> std::pair<std::vector<double>, std::vector<double>> {
            Vector1d1 mx, mn;
            HGP_3D_Mesh_Curvature_C2(to_vec3d1(vecs), face_ids, mx, mn);
            return {mx, mn};
        }, py::arg("vecs"), py::arg("face_ids"),
        "网格曲率C2（face_ids为Nx3整数列表）。返回 (max_curs, min_curs)。");

    // HGP_3D_Mesh_Curvature_C3（out: max/min_curs + directions）
    m.def("HGP_3D_Mesh_Curvature_C3",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2)
           -> std::tuple<std::vector<double>, std::vector<double>,
                         std::vector<std::vector<double>>, std::vector<std::vector<double>>> {
            Vector1d1 mx, mn; Vector3d1 dmx, dmn;
            HGP_3D_Mesh_Curvature_C3(to_vec3d1(vecs), f0, f1, f2, mx, mn, dmx, dmn);
            return {mx, mn, from_vec3d1(dmx), from_vec3d1(dmn)};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "网格曲率C3。返回 (max_curs, min_curs, max方向, min方向)。");

    m.def("HGP_3D_Mesh_Curvature_C4",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<std::vector<int>>& face_ids)
           -> std::tuple<std::vector<double>, std::vector<double>,
                         std::vector<std::vector<double>>, std::vector<std::vector<double>>> {
            Vector1d1 mx, mn; Vector3d1 dmx, dmn;
            HGP_3D_Mesh_Curvature_C4(to_vec3d1(vecs), face_ids, mx, mn, dmx, dmn);
            return {mx, mn, from_vec3d1(dmx), from_vec3d1(dmn)};
        }, py::arg("vecs"), py::arg("face_ids"),
        "网格曲率C4（face_ids Nx3）。返回 (max, min, max方向, min方向)。");

    m.def("HGP_3D_Mesh_Curvature_C5",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2)
           -> std::tuple<std::vector<double>, std::vector<double>,
                         std::vector<std::vector<double>>, std::vector<std::vector<double>>,
                         std::vector<std::vector<double>>> {
            Vector1d1 mx, mn; Vector3d1 dmx, dmn, normals;
            HGP_3D_Mesh_Curvature_C5(to_vec3d1(vecs), f0, f1, f2, mx, mn, dmx, dmn, normals);
            return {mx, mn, from_vec3d1(dmx), from_vec3d1(dmn), from_vec3d1(normals)};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "网格曲率C5（含法向量）。");

    m.def("HGP_3D_Mesh_Curvature_C6",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<std::vector<int>>& face_ids)
           -> std::tuple<std::vector<double>, std::vector<double>,
                         std::vector<std::vector<double>>, std::vector<std::vector<double>>,
                         std::vector<std::vector<double>>> {
            Vector1d1 mx, mn; Vector3d1 dmx, dmn, normals;
            HGP_3D_Mesh_Curvature_C6(to_vec3d1(vecs), face_ids, mx, mn, dmx, dmn, normals);
            return {mx, mn, from_vec3d1(dmx), from_vec3d1(dmn), from_vec3d1(normals)};
        }, py::arg("vecs"), py::arg("face_ids"),
        "网格曲率C6（face_ids Nx3，含法向量）。");

    // HGP_3D_Triangle_Mesh_Boundary_C1 / C2（out: bools）
    m.def("HGP_3D_Triangle_Mesh_Boundary_C1",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2)
           -> std::vector<bool> {
            Vector1b1 bools;
            HGP_3D_Triangle_Mesh_Boundary_C1(to_vec3d1(vecs), f0, f1, f2, bools);
            return bools;
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "标记三角网格边界顶点（C1）。返回 bool 列表（每个顶点是否在边界）。");

    m.def("HGP_3D_Triangle_Mesh_Boundary_C2",
        [](const std::string& path) -> std::vector<bool> {
            Vector1b1 bools;
            HGP_3D_Triangle_Mesh_Boundary_C2(path.c_str(), bools);
            return bools;
        }, py::arg("path"), "标记网格边界顶点（C2，文件）。");

    // HGP_3D_Connecting_Segments_C1（inout: segments -> lines，2D）
    m.def("HGP_3D_Connecting_Segments_C1",
        [](std::vector<std::vector<std::vector<double>>> segs)
           -> std::pair<std::vector<std::vector<std::vector<double>>>,
                        std::vector<std::vector<std::vector<double>>>> {
            Vector2d2 s = to_vec2d2(segs), lines;
            HGP_3D_Connecting_Segments_C1(s, lines);
            return {from_vec2d2(s), from_vec2d2(lines)};
        }, py::arg("segments"),
        "2D线段连接（C1）。返回 (更新后的segments, 连接后的折线列表)。");

    // HGP_3D_Connecting_Segments_C2（3D）
    m.def("HGP_3D_Connecting_Segments_C2",
        [](std::vector<std::vector<std::vector<double>>> segs)
           -> std::pair<std::vector<std::vector<std::vector<double>>>,
                        std::vector<std::vector<std::vector<double>>>> {
            Vector3d2 s = to_vec3d2(segs), lines;
            HGP_3D_Connecting_Segments_C2(s, lines);
            return {from_vec3d2(s), from_vec3d2(lines)};
        }, py::arg("segments"),
        "3D线段连接（C2）。返回 (更新后的segments, 连接后的折线列表)。");

    // HGP_3D_Triangle_Mesh_Boundary_C3/C4/C5（out: boundaries Vector3d2）
    m.def("HGP_3D_Triangle_Mesh_Boundary_C3",
        [](std::vector<std::vector<double>> vecs,
           std::vector<int> f0, std::vector<int> f1, std::vector<int> f2)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>,
                         std::vector<std::vector<std::vector<double>>>> {
            Vector3d1 v = to_vec3d1(vecs); Vector3d2 bounds;
            HGP_3D_Triangle_Mesh_Boundary_C3(v, f0, f1, f2, bounds);
            return {from_vec3d1(v), f0, f1, f2, from_vec3d2(bounds)};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "提取网格边界折线C3（inout vecs/faces）。");

    m.def("HGP_3D_Triangle_Mesh_Boundary_C4",
        [](std::vector<std::vector<double>> vecs,
           std::vector<std::vector<int>> face_ids)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<std::vector<int>>,
                         std::vector<std::vector<std::vector<double>>>> {
            Vector3d1 v = to_vec3d1(vecs); Vector3d2 bounds;
            HGP_3D_Triangle_Mesh_Boundary_C4(v, face_ids, bounds);
            return {from_vec3d1(v), face_ids, from_vec3d2(bounds)};
        }, py::arg("vecs"), py::arg("face_ids"),
        "提取网格边界折线C4（face_ids Nx3）。");

    m.def("HGP_3D_Triangle_Mesh_Boundary_C5",
        [](const std::string& path)
           -> std::vector<std::vector<std::vector<double>>> {
            Vector3d2 bounds;
            HGP_3D_Triangle_Mesh_Boundary_C5(path.c_str(), bounds);
            return from_vec3d2(bounds);
        }, py::arg("path"), "提取网格边界折线C5（文件）。返回边界折线列表。");

    // HGP_Mesh_Laplace_Smooth_C1（文件IO）
    m.def("HGP_Mesh_Laplace_Smooth_C1",
        [](const std::string& in, const std::string& out, int nb){
            HGP_Mesh_Laplace_Smooth_C1(in.c_str(), out.c_str(), nb);
        }, py::arg("in_path"), py::arg("out_path"), py::arg("laplace_nb"),
        "Laplace 平滑网格（C1，文件IO）。");

    // HGP_3D_Triangle_Mesh_Vecs_Neighbors（out: neighs Vector1i2）
    m.def("HGP_3D_Triangle_Mesh_Vecs_Neighbors",
        [](std::vector<std::vector<double>> vecs,
           std::vector<int> f0, std::vector<int> f1, std::vector<int> f2)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>,
                         std::vector<std::vector<int>>> {
            Vector3d1 v = to_vec3d1(vecs); Vector1i2 neighs;
            HGP_3D_Triangle_Mesh_Vecs_Neighbors(v, f0, f1, f2, neighs);
            return {from_vec3d1(v), f0, f1, f2, neighs};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "计算顶点邻居。返回 (vecs, f0, f1, f2, neighs[顶点->邻居顶点列表])。");

    // HGP_Mesh_Laplace_Smooth_C2（inout vecs, faces）
    m.def("HGP_Mesh_Laplace_Smooth_C2",
        [](std::vector<std::vector<double>> vecs,
           std::vector<int> f0, std::vector<int> f1, std::vector<int> f2, int nb)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>> {
            Vector3d1 v = to_vec3d1(vecs);
            HGP_Mesh_Laplace_Smooth_C2(v, f0, f1, f2, nb);
            return {from_vec3d1(v), f0, f1, f2};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
           py::arg("laplace_nb"),
        "Laplace 平滑（C2，原地修改顶点）。");

    // HGP_3D_Triangle_Mesh_Vecs_Faces（out: surface_vectices_to_face Vector1i2）
    m.def("HGP_3D_Triangle_Mesh_Vecs_Faces",
        [](std::vector<std::vector<double>> vecs,
           std::vector<int> f0, std::vector<int> f1, std::vector<int> f2)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>,
                         std::vector<std::vector<int>>> {
            Vector3d1 v = to_vec3d1(vecs); Vector1i2 vtf;
            HGP_3D_Triangle_Mesh_Vecs_Faces(v, f0, f1, f2, vtf);
            return {from_vec3d1(v), f0, f1, f2, vtf};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "计算顶点-面映射。返回每个顶点所属的面索引列表。");

    // HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges（out: Vector1i3）
    m.def("HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges",
        [](std::vector<std::vector<double>> vecs,
           std::vector<int> f0, std::vector<int> f1, std::vector<int> f2)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>,
                         std::vector<std::vector<std::vector<int>>>> {
            Vector3d1 v = to_vec3d1(vecs); Vector1i3 vne;
            HGP_3D_Triangle_Mesh_Vecs_Neighbor_Edges(v, f0, f1, f2, vne);
            return {from_vec3d1(v), f0, f1, f2, vne};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "计算顶点邻居边映射（3层嵌套整数列表）。");

    // HGP_Mesh_Laplace_Smooth_by_Curvature（inout vecs/faces + out low_curvature）
    m.def("HGP_Mesh_Laplace_Smooth_by_Curvature",
        [](std::vector<std::vector<double>> vecs,
           std::vector<int> f0, std::vector<int> f1, std::vector<int> f2)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>, double> {
            Vector3d1 v = to_vec3d1(vecs); double lc=0.0;
            HGP_Mesh_Laplace_Smooth_by_Curvature(v, f0, f1, f2, lc);
            return {from_vec3d1(v), f0, f1, f2, lc};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "基于曲率的 Laplace 平滑。返回 (新顶点, f0, f1, f2, low_curvature)。");

    // HGP_Mesh_Loop_Subdivision_Own_Version（文件IO）
    m.def("HGP_Mesh_Loop_Subdivision_Own_Version",
        [](const std::string& in, int step, const std::string& out, int nb){
            HGP_Mesh_Loop_Subdivision_Own_Version(in.c_str(), step, out.c_str(), nb);
        }, py::arg("in_path"), py::arg("step"), py::arg("out_path"), py::arg("laplace_nb"),
        "自定义 Loop 细分（含 Laplace 平滑，文件IO）。");

    // HGP_Rotation_Obj（文件IO，axis 为 Vector3d）
    m.def("HGP_Rotation_Obj",
        [](const std::string& path, double angle, const std::vector<double>& axis){
            HGP_Rotation_Obj(path.c_str(), angle, to_vec3d(axis));
        }, py::arg("path"), py::arg("angle"), py::arg("axis"),
        "旋转 .obj 网格文件（原地修改）。angle 单位：弧度。");

    // HGP_Slicer_Mesh（out: offsetses Vector3d3, offsets Vector3d2）
    m.def("HGP_Slicer_Mesh",
        [](const std::string& path,
           const std::vector<double>& plane_normal,
           const std::vector<double>& plane_d)
           -> std::pair<std::vector<std::vector<std::vector<std::vector<double>>>>,
                        std::vector<std::vector<std::vector<double>>>> {
            Vector3d3 offsetses; Vector3d2 offsets;
            HGP_Slicer_Mesh(path.c_str(), to_vec3d(plane_normal), plane_d,
                            offsetses, offsets);
            return {from_vec3d3(offsetses), from_vec3d2(offsets)};
        }, py::arg("path"), py::arg("plane_normal"), py::arg("plane_d"),
        "网格切片。返回 (offsetses[每层多条轮廓], offsets[合并轮廓])。");

    // HGP_Shortest_Geodesic_Path_C1（out: xyzs）
    m.def("HGP_Shortest_Geodesic_Path_C1",
        [](const std::string& path) -> std::vector<std::vector<double>> {
            Vector3d1 xyzs;
            HGP_Shortest_Geodesic_Path_C1(path.c_str(), xyzs);
            return from_vec3d1(xyzs);
        }, py::arg("path"), "计算网格上最短测地路径（C1）。");

    // HGP_Shortest_Geodesic_Path_C3（source, target 按值传入）
    m.def("HGP_Shortest_Geodesic_Path_C3",
        [](const std::string& path,
           const std::vector<double>& source,
           const std::vector<double>& target) -> std::vector<std::vector<double>> {
            Vector3d1 output;
            HGP_Shortest_Geodesic_Path_C3(
                path.c_str(), to_vec3d(source), to_vec3d(target), output);
            return from_vec3d1(output);
        }, py::arg("path"), py::arg("source"), py::arg("target"),
        "计算网格上两点间最短测地路径（C3）。");

    // HGP_Shortest_Geodesic_Path_C4（sources, targets 按值传入，out: xyzes Vector3d2）
    m.def("HGP_Shortest_Geodesic_Path_C4",
        [](const std::string& path,
           const std::vector<std::vector<double>>& sources,
           const std::vector<std::vector<double>>& targets)
           -> std::vector<std::vector<std::vector<double>>> {
            Vector3d2 xyzes;
            HGP_Shortest_Geodesic_Path_C4(
                path.c_str(), to_vec3d1(sources), to_vec3d1(targets), xyzes);
            return from_vec3d2(xyzes);
        }, py::arg("path"), py::arg("sources"), py::arg("targets"),
        "批量最短测地路径（C4）。返回每对 (source, target) 的路径点列表。");

    // HGP_Geodesic_Distance
    m.def("HGP_Geodesic_Distance",
        [](const std::string& path,
           const std::vector<double>& source,
           const std::vector<double>& target) -> double {
            return HGP_Geodesic_Distance(
                path.c_str(), to_vec3d(source), to_vec3d(target));
        }, py::arg("path"), py::arg("source"), py::arg("target"),
        "计算网格上两点间测地距离。");

    // HGP_Project_Points_Onto_Surface_C1 / C2
    m.def("HGP_Project_Points_Onto_Surface_C1",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2,
           const std::vector<std::vector<double>>& points)
           -> std::vector<std::vector<double>> {
            return from_vec3d1(HGP_Project_Points_Onto_Surface_C1(
                to_vec3d1(vecs), f0, f1, f2, to_vec3d1(points)));
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
           py::arg("points"),
        "将点集投影到三角网格表面（C1）。");

    m.def("HGP_Project_Points_Onto_Surface_C2",
        [](const std::string& path,
           const std::vector<std::vector<double>>& points)
           -> std::vector<std::vector<double>> {
            return from_vec3d1(HGP_Project_Points_Onto_Surface_C2(
                path.c_str(), to_vec3d1(points)));
        }, py::arg("path"), py::arg("points"),
        "将点集投影到网格表面（C2，文件）。");

    // HGP_3D_Triangel_Mesh_Most_Inside_Point（out: inside Vector3d）
    m.def("HGP_3D_Triangel_Mesh_Most_Inside_Point",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2)
           -> std::vector<double> {
            Vector3d inside;
            HGP_3D_Triangel_Mesh_Most_Inside_Point(to_vec3d1(vecs), f0, f1, f2, inside);
            return from_vec3d(inside);
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "返回三角网格的最深内部点 [x,y,z]。");

    // HGP_3D_Triangle_Mesh_Area
    m.def("HGP_3D_Triangle_Mesh_Area",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2)
           -> double {
            return HGP_3D_Triangle_Mesh_Area(to_vec3d1(vecs), f0, f1, f2);
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "计算三角网格总面积。");

    // HGP_3D_Convex_Hulls_C1 / C2 / C3 / C4
    m.def("HGP_3D_Convex_Hulls_C1",
        [](const std::vector<std::vector<double>>& vec)
           -> std::vector<std::vector<double>> {
            Vector3d1 hull;
            HGP_3D_Convex_Hulls_C1(to_vec3d1(vec), hull);
            return from_vec3d1(hull);
        }, py::arg("vec"), "3D凸包顶点（C1）。");

    m.def("HGP_3D_Convex_Hulls_C2",
        [](const std::vector<std::vector<double>>& vec)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>> {
            Vector3d1 hull; Vector1i1 s0, s1, s2;
            HGP_3D_Convex_Hulls_C2(to_vec3d1(vec), hull, s0, s1, s2);
            return {from_vec3d1(hull), s0, s1, s2};
        }, py::arg("vec"),
        "3D凸包（C2）。返回 (凸包顶点, 面索引0/1/2)。");

    m.def("HGP_3D_Convex_Hulls_C3",
        [](const std::vector<std::vector<double>>& vec)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<std::vector<double>>,
                         std::vector<std::vector<double>>> {
            Vector3d1 hull, plane_p, plane_n;
            HGP_3D_Convex_Hulls_C3(to_vec3d1(vec), hull, plane_p, plane_n);
            return {from_vec3d1(hull), from_vec3d1(plane_p), from_vec3d1(plane_n)};
        }, py::arg("vec"),
        "3D凸包（C3）。返回 (凸包顶点, 面平面点列表, 面法向量列表)。");

    m.def("HGP_3D_Convex_Hulls_C4",
        [](const std::vector<std::vector<double>>& vec)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>, std::vector<int>, std::vector<int>,
                         std::vector<std::vector<double>>,
                         std::vector<std::vector<double>>> {
            Vector3d1 hull, pp, pn; Vector1i1 s0, s1, s2;
            HGP_3D_Convex_Hulls_C4(to_vec3d1(vec), hull, s0, s1, s2, pp, pn);
            return {from_vec3d1(hull), s0, s1, s2, from_vec3d1(pp), from_vec3d1(pn)};
        }, py::arg("vec"),
        "3D凸包（C4）。返回 (凸包顶点, 面索引0/1/2, 平面点, 平面法向)。");

    // HGP_Mesh_Field_Query_C1 / C2 / C3
    m.def("HGP_Mesh_Field_Query_C1",
        [](const std::string& path,
           const std::vector<std::vector<double>>& gradients,
           const std::vector<std::vector<double>>& input_pts)
           -> std::vector<std::vector<double>> {
            Vector3d1 out;
            HGP_Mesh_Field_Query_C1(
                path.c_str(), to_vec3d1(gradients), to_vec3d1(input_pts), out);
            return from_vec3d1(out);
        }, py::arg("path"), py::arg("gradients"), py::arg("input_points"),
        "网格场查询C1（梯度向量插值）。");

    m.def("HGP_Mesh_Field_Query_C2",
        [](const std::string& path,
           const std::vector<double>& gradient_values,
           const std::vector<std::vector<double>>& input_pts) -> std::vector<double> {
            Vector1d1 out;
            HGP_Mesh_Field_Query_C2(
                path.c_str(), gradient_values, to_vec3d1(input_pts), out);
            return out;
        }, py::arg("path"), py::arg("gradient_values"), py::arg("input_points"),
        "网格场查询C2（标量插值）。");

    m.def("HGP_Mesh_Field_Query_C3",
        [](const std::string& path,
           const std::vector<double>& gradient_values,
           const std::vector<std::vector<std::vector<double>>>& input_pt_es)
           -> std::vector<std::vector<double>> {
            Vector1d2 out;
            HGP_Mesh_Field_Query_C3(
                path.c_str(), gradient_values, to_vec3d2(input_pt_es), out);
            return out;
        }, py::arg("path"), py::arg("gradient_values"), py::arg("input_point_es"),
        "网格场查询C3（多组标量插值）。");

    // HGP_Curvature_Mesh（out: max/min curs + directions）
    m.def("HGP_Curvature_Mesh",
        [](const std::string& path,
           const std::vector<std::vector<double>>& input_pts)
           -> std::tuple<std::vector<double>, std::vector<double>,
                         std::vector<std::vector<double>>, std::vector<std::vector<double>>> {
            Vector1d1 mx, mn; Vector3d1 dmx, dmn;
            HGP_Curvature_Mesh(path.c_str(), to_vec3d1(input_pts), mx, mn, dmx, dmn);
            return {mx, mn, from_vec3d1(dmx), from_vec3d1(dmn)};
        }, py::arg("path"), py::arg("input_points"),
        "网格文件曲率查询。返回 (max_curs, min_curs, max方向, min方向)。");

    // HGP_Normal_Mesh_C1 / C2
    m.def("HGP_Normal_Mesh_C1",
        [](const std::string& path,
           const std::vector<std::vector<double>>& pts)
           -> std::vector<std::vector<double>> {
            Vector3d1 normals;
            HGP_Normal_Mesh_C1(path.c_str(), to_vec3d1(pts), normals);
            return from_vec3d1(normals);
        }, py::arg("path"), py::arg("mesh_points"),
        "查询网格上指定点的法向量（C1）。");

    m.def("HGP_Normal_Mesh_C2",
        [](const std::string& path,
           const std::vector<std::vector<std::vector<double>>>& ptses)
           -> std::vector<std::vector<std::vector<double>>> {
            Vector3d2 normalses;
            HGP_Normal_Mesh_C2(path.c_str(), to_vec3d2(ptses), normalses);
            return from_vec3d2(normalses);
        }, py::arg("path"), py::arg("mesh_pointses"),
        "批量查询网格上点的法向量（C2，多组）。");

    // HGP_3D_Mesh_Normal_C1 / C2（out: normals）
    m.def("HGP_3D_Mesh_Normal_C1",
        [](const std::vector<std::vector<double>>& ps,
           const std::vector<std::vector<int>>& face_ids)
           -> std::vector<std::vector<double>> {
            Vector3d1 normals;
            HGP_3D_Mesh_Normal_C1(to_vec3d1(ps), face_ids, normals);
            return from_vec3d1(normals);
        }, py::arg("ps"), py::arg("face_ids"),
        "计算网格顶点法向量（C1，face_ids Nx3）。");

    m.def("HGP_3D_Mesh_Normal_C2",
        [](const std::vector<std::vector<double>>& ps,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2)
           -> std::vector<std::vector<double>> {
            Vector3d1 normals;
            HGP_3D_Mesh_Normal_C2(to_vec3d1(ps), f0, f1, f2, normals);
            return from_vec3d1(normals);
        }, py::arg("ps"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "计算网格顶点法向量（C2，分量传参）。");

    // HGP_3D_Mesh_Center_C1 / C2
    m.def("HGP_3D_Mesh_Center_C1",
        [](const std::vector<std::vector<std::vector<double>>>& ps) -> std::vector<double> {
            return from_vec3d(HGP_3D_Mesh_Center_C1(to_vec3d2(ps)));
        }, py::arg("ps"), "多组3D点的重心（C1，Vector3d2）。");

    m.def("HGP_3D_Mesh_Center_C2",
        [](const std::vector<std::vector<double>>& ps) -> std::vector<double> {
            return from_vec3d(HGP_3D_Mesh_Center_C2(to_vec3d1(ps)));
        }, py::arg("ps"), "3D点集的重心（C2，Vector3d1）。");

    // HGP_3D_Mesh_Boundingbox_C1 / C2（out: min/max corner）
    m.def("HGP_3D_Mesh_Boundingbox_C1",
        [](const std::vector<std::vector<std::vector<double>>>& ps)
           -> std::pair<std::vector<double>, std::vector<double>> {
            Vector3d mn, mx;
            HGP_3D_Mesh_Boundingbox_C1(to_vec3d2(ps), mn, mx);
            return {from_vec3d(mn), from_vec3d(mx)};
        }, py::arg("ps"), "多组3D点的包围盒（C1）。返回 (min, max)。");

    m.def("HGP_3D_Mesh_Boundingbox_C2",
        [](const std::vector<std::vector<double>>& ps)
           -> std::pair<std::vector<double>, std::vector<double>> {
            Vector3d mn, mx;
            HGP_3D_Mesh_Boundingbox_C2(to_vec3d1(ps), mn, mx);
            return {from_vec3d(mn), from_vec3d(mx)};
        }, py::arg("ps"), "3D点集包围盒（C2）。返回 (min, max)。");

    // HGP_Surface_Decomposition（out: face_sdf, regions_nb, face_segments）
    m.def("HGP_Surface_Decomposition",
        [](const std::string& path)
           -> std::tuple<std::vector<double>, int, std::vector<int>> {
            Vector1d1 sdf; int nb=0; Vector1i1 segs;
            HGP_Surface_Decomposition(path.c_str(), sdf, nb, segs);
            return {sdf, nb, segs};
        }, py::arg("path"),
        "网格表面分解。返回 (face_sdf, regions_nb, face_segments)。");

    // HGP_3D_Mesh_Gradient（out: vecs_gradients, faces_gradients）
    m.def("HGP_3D_Mesh_Gradient",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2,
           const std::vector<double>& psd)
           -> std::pair<std::vector<std::vector<double>>,
                        std::vector<std::vector<double>>> {
            Vector3d1 vg, fg;
            HGP_3D_Mesh_Gradient(to_vec3d1(vecs), f0, f1, f2, psd, vg, fg);
            return {from_vec3d1(vg), from_vec3d1(fg)};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
           py::arg("psd"),
        "计算网格梯度。返回 (顶点梯度, 面梯度)，各为 [[x,y,z], ...]。");

    // HGP_Intergral_Curvature（out: output_points, output_rates）
    m.def("HGP_Intergral_Curvature",
        [](const std::vector<std::vector<double>>& input_pts,
           int nb, double radius, double thresholder)
           -> std::pair<std::vector<std::vector<double>>, std::vector<double>> {
            Vector2d1 out_pts; Vector1d1 out_rates;
            HGP_Intergral_Curvature(
                to_vec2d1(input_pts), nb, radius, thresholder, out_pts, out_rates);
            return {from_vec2d1(out_pts), out_rates};
        }, py::arg("input_points"), py::arg("sampling_points_nb"),
           py::arg("radius"), py::arg("thresholder"),
        "2D��分曲率。返回 (output_points, output_rates)。");

    // HGP_3D_Mesh_Extract_Isoline（out: isolines Vector3d2）
    m.def("HGP_3D_Mesh_Extract_Isoline",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0, const std::vector<int>& f1, const std::vector<int>& f2,
           const std::vector<double>& psd, double d)
           -> std::pair<bool, std::vector<std::vector<std::vector<double>>>> {
            Vector3d2 isolines;
            bool ok = HGP_3D_Mesh_Extract_Isoline(
                to_vec3d1(vecs), f0, f1, f2, psd, d, isolines);
            return {ok, from_vec3d2(isolines)};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
           py::arg("psd"), py::arg("d"),
        "提取网格等值线。返回 (bool, 等值线列表[条][点][xyz])。");

    // HGP_BSplineCurveFit（out: output）
    m.def("HGP_BSplineCurveFit",
        [](const std::vector<std::vector<double>>& samples)
           -> std::vector<std::vector<double>> {
            Vector3d1 output;
            HGP_BSplineCurveFit(to_vec3d1(samples), output);
            return from_vec3d1(output);
        }, py::arg("samples"), "3D点集 B 样条曲线拟合，返回拟合后的曲线点列表。");

    // HGP_Cut_Surface / HGP_Cut_Surface_by_Multi_Boundaries（output_path 为 char* out-buf）
    m.def("HGP_Cut_Surface",
        [](const std::vector<std::vector<double>>& boundary,
           const std::vector<double>& inside_point,
           const std::string& full_path) -> std::string {
            char buf[1024] = {};
            HGP_Cut_Surface(to_vec3d1(boundary), to_vec3d(inside_point),
                            full_path.c_str(), buf);
            return std::string(buf);
        }, py::arg("boundary"), py::arg("inside_point"), py::arg("full_path"),
        "按单条边界裁剪网格。返回输出文件路径字符串。");

        // HGP_Cut_Surface_by_Multi_Boundaries
    //   multi_boundary: list[list[[x,y,z]]]  —— 多条边界曲线
    //   inside_point : [x,y,z]
    //   full_path    : 输入网格文件路径
    //   返回         : 输出文件路径字符串
    m.def("HGP_Cut_Surface_by_Multi_Boundaries",
        [](const std::vector<std::vector<std::vector<double>>>& multi_boundary,
           const std::vector<double>& inside_point,
           const std::string& full_path) -> std::string {
            char buf[1024] = {};
            HGP_Cut_Surface_by_Multi_Boundaries(
                to_vec3d2(multi_boundary),
                to_vec3d(inside_point),
                full_path.c_str(),
                buf);
            return std::string(buf);
        }, py::arg("multi_boundary"), py::arg("inside_point"), py::arg("full_path"),
        "按多条边界裁剪网格。返回输出文件路径字符串。");
        // ════════════════════════════════════════════════════
   
    // ── 补全：缺失的 HGP 函数
    // ════════════════════════════════════════════════════
    // ── HGP2D: Distance ─────────────────────────────────

    m.def("HGP_2D_Distance_Point_Point",
        [](const std::vector<double>& p0, const std::vector<double>& p1) -> double {
            return HGP_2D_Distance_Point_Point(to_vec2d(p0), to_vec2d(p1));
        }, py::arg("p_0"), py::arg("p_1"), "2D两点距离。");

    m.def("HGP_2D_Distance_Point_Line",
        [](const std::vector<double>& v,
           const std::vector<double>& l0, const std::vector<double>& l1) -> double {
            return HGP_2D_Distance_Point_Line(to_vec2d(v), to_vec2d(l0), to_vec2d(l1));
        }, py::arg("v"), py::arg("l_0"), py::arg("l_1"), "2D点到直线距离。");

    m.def("HGP_2D_Distance_Point_Segment",
        [](const std::vector<double>& v,
           const std::vector<double>& s0, const std::vector<double>& s1) -> double {
            return HGP_2D_Distance_Point_Segment(to_vec2d(v), to_vec2d(s0), to_vec2d(s1));
        }, py::arg("v"), py::arg("s_0"), py::arg("s_1"), "2D点到线段距离。");

    m.def("HGP_2D_Distance_Segment_Segment",
        [](const std::vector<double>& s0, const std::vector<double>& s1,
           const std::vector<double>& e0, const std::vector<double>& e1) -> double {
            return HGP_2D_Distance_Segment_Segment(
                to_vec2d(s0), to_vec2d(s1), to_vec2d(e0), to_vec2d(e1));
        }, py::arg("s_0"), py::arg("s_1"), py::arg("e_0"), py::arg("e_1"),
        "2D两线段距离。");

    m.def("HGP_2D_Distance_Point_Polygon",
        [](const std::vector<double>& p,
           const std::vector<std::vector<double>>& py_) -> double {
            return HGP_2D_Distance_Point_Polygon(to_vec2d(p), to_vec2d1(py_));
        }, py::arg("p"), py::arg("py"), "2D点到多边形距离。");

    // ── HGP2D: Location ─────────────────────────────────

    m.def("HGP_2D_Location_Point_Polygon",
        [](const std::vector<double>& p,
           const std::vector<std::vector<double>>& py_) -> bool {
            return HGP_2D_Location_Point_Polygon(to_vec2d(p), to_vec2d1(py_));
        }, py::arg("p"), py::arg("py"), "判断2D点是否在多边形内。");

    m.def("HGP_2D_Location_Points_Polygon",
        [](const std::vector<std::vector<double>>& ps,
           const std::vector<std::vector<double>>& py_) -> bool {
            return HGP_2D_Location_Points_Polygon(to_vec2d1(ps), to_vec2d1(py_));
        }, py::arg("ps"), py::arg("py"), "判断2D点集是否均在多边形内。");

    // ── HGP2D: Sampling ─────────────────────────────────

    m.def("HGP_2D_Polygon_Dart_Sampling",
        [](const std::vector<std::vector<double>>& py_, double d, int total_iter)
           -> std::vector<std::vector<double>> {
            Vector2d1 sampling_points;
            HGP_2D_Polygon_Dart_Sampling(to_vec2d1(py_), d, sampling_points, total_iter);
            return from_vec2d1(sampling_points);
        }, py::arg("py"), py::arg("d"), py::arg("total_iter"),
        "多边形 dart 采样。d 为包围盒对角线百分比。返回采样点列表。");

    m.def("HGP_2D_Polygon_Regular_Sampling_C1",
        [](const std::vector<std::vector<double>>& py_, double d)
           -> std::vector<std::vector<double>> {
            return from_vec2d1(HGP_2D_Polygon_Regular_Sampling_C1(to_vec2d1(py_), d));
        }, py::arg("py"), py::arg("d"), "多边形规则采样C1。");

    m.def("HGP_2D_Polygon_Regular_Sampling_C2",
        [](const std::vector<std::vector<double>>& py_, double d)
           -> std::pair<std::vector<std::vector<double>>,
                        std::vector<std::pair<int,int>>> {
            VectorPI1 nb;
            auto pts = HGP_2D_Polygon_Regular_Sampling_C2(to_vec2d1(py_), d, nb);
            return {from_vec2d1(pts), from_pi1(nb)};
        }, py::arg("py"), py::arg("d"),
        "多边形规则采样C2。返回 (采样点, 邻居对列表)。");

    m.def("HGP_2D_Polygon_Regular_Sampling_C3",
        [](const std::vector<std::vector<double>>& py_, double d, bool compute_neighbors)
           -> std::pair<std::vector<std::vector<double>>,
                        std::vector<std::pair<int,int>>> {
            VectorPI1 nb;
            auto pts = HGP_2D_Polygon_Regular_Sampling_C3(
                to_vec2d1(py_), d, nb, compute_neighbors);
            return {from_vec2d1(pts), from_pi1(nb)};
        }, py::arg("py"), py::arg("d"), py::arg("compute_neighbors"),
        "多边形规则采样C3。返回 (采样点, 邻居对列表)。");

    // ── HGP2D: Intersection ─────────────────────────────

    m.def("HGP_2D_Intersection_Segment_Segment",
        [](const std::vector<double>& s0s, const std::vector<double>& s0e,
           const std::vector<double>& s1s, const std::vector<double>& s1e)
           -> std::pair<bool, std::vector<double>> {
            Vector2d inter;
            bool hit = HGP_2D_Intersection_Segment_Segment(
                to_vec2d(s0s), to_vec2d(s0e),
                to_vec2d(s1s), to_vec2d(s1e), inter);
            return {hit, from_vec2d(inter)};
        }, py::arg("s_0_s"), py::arg("s_0_e"), py::arg("s_1_s"), py::arg("s_1_e"),
        "2D两线段求交。返回 (是否相交, 交点[x,y])。");

    m.def("HGP_2D_Intersection_Line_Line",
        [](const std::vector<double>& s0s, const std::vector<double>& s0e,
           const std::vector<double>& s1s, const std::vector<double>& s1e)
           -> std::pair<bool, std::vector<double>> {
            Vector2d inter;
            bool hit = HGP_2D_Intersection_Line_Line(
                to_vec2d(s0s), to_vec2d(s0e),
                to_vec2d(s1s), to_vec2d(s1e), inter);
            return {hit, from_vec2d(inter)};
        }, py::arg("s_0_s"), py::arg("s_0_e"), py::arg("s_1_s"), py::arg("s_1_e"),
        "2D两直线求交。返回 (是否相交, 交点[x,y])。");

    // ── HGP2D: Polygon properties ───────────────────────

    m.def("HGP_2D_Polygon_Is_Clockwise_Oriented",
        [](const std::vector<std::vector<double>>& ps) -> bool {
            return HGP_2D_Polygon_Is_Clockwise_Oriented(to_vec2d1(ps));
        }, py::arg("ps"), "判断2D多边形是否为顺时针方向。");

    m.def("HGP_2D_Polygon_Area",
        [](const std::vector<std::vector<double>>& py_) -> double {
            return HGP_2D_Polygon_Area(to_vec2d1(py_));
        }, py::arg("py"), "2D多边形面积。");

    m.def("HGP_2D_Polygon_One_Offsets",
        [](const std::vector<std::vector<double>>& poly, double d)
           -> std::vector<std::vector<std::vector<double>>> {
            Vector2d2 offset_polys;
            HGP_2D_Polygon_One_Offsets(to_vec2d1(poly), d, offset_polys);
            return from_vec2d2(offset_polys);
        }, py::arg("poly"), py::arg("d"),
        "单个多边形偏置（offset）。返回偏置后多边形列表。");

    m.def("HGP_2D_Nearest_Point_Polygon_C1",
        [](const std::vector<double>& v,
           const std::vector<std::vector<double>>& poly) -> std::vector<double> {
            return from_vec2d(HGP_2D_Nearest_Point_Polygon_C1(to_vec2d(v), to_vec2d1(poly)));
        }, py::arg("v"), py::arg("poly"),
        "多边形上距 v 最近的点（C1）。返回 [x,y]。");

    m.def("HGP_2d_Polygon_Boundingbox",
        [](const std::vector<std::vector<double>>& ps)
           -> std::pair<std::vector<double>, std::vector<double>> {
            Vector2d mn, mx;
            HGP_2d_Polygon_Boundingbox(to_vec2d1(ps), mn, mx);
            return {from_vec2d(mn), from_vec2d(mx)};
        }, py::arg("ps"), "2D多边形包围盒。返回 (min_corner[x,y], max_corner[x,y])。");

    m.def("HGP_2D_Convex_Hulls",
        [](const std::vector<std::vector<double>>& vec)
           -> std::vector<std::vector<double>> {
            Vector2d1 hull;
            HGP_2D_Convex_Hulls(to_vec2d1(vec), hull);
            return from_vec2d1(hull);
        }, py::arg("vec"), "2D点集凸包。返回凸包顶点列表。");

    // ── HGP3D: 补全缺失函数 ──────────────────────────────

    m.def("HGP_3D_Distance_Point_Point",
        [](const std::vector<double>& v0, const std::vector<double>& v1) -> double {
            return HGP_3D_Distance_Point_Point(to_vec3d(v0), to_vec3d(v1));
        }, py::arg("v_0"), py::arg("v_1"), "3D两点距离。");

    m.def("HGP_3D_Distance_Point_Segment",
        [](const std::vector<double>& p,
           const std::vector<double>& ss, const std::vector<double>& se) -> double {
            return HGP_3D_Distance_Point_Segment(to_vec3d(p), to_vec3d(ss), to_vec3d(se));
        }, py::arg("p"), py::arg("s_s"), py::arg("s_e"), "3D点到线段距离。");

    m.def("HGP_3D_Plane_Fitting",
        [](const std::vector<std::vector<double>>& points)
           -> std::pair<std::vector<double>, std::vector<double>> {
            Vector3d plane_p, plane_n;
            HGP_3D_Plane_Fitting(to_vec3d1(points), plane_p, plane_n);
            return {from_vec3d(plane_p), from_vec3d(plane_n)};
        }, py::arg("points"),
        "3D点集平面拟合。返回 (平面点[x,y,z], 平面法向[x,y,z])。");

    m.def("HGP_3D_Intersection_Segment_Plane",
        [](const std::vector<double>& ss, const std::vector<double>& se,
           const std::vector<double>& pp, const std::vector<double>& pn)
           -> std::pair<bool, std::vector<double>> {
            Vector3d inter;
            bool hit = HGP_3D_Intersection_Segment_Plane(
                to_vec3d(ss), to_vec3d(se), to_vec3d(pp), to_vec3d(pn), inter);
            return {hit, from_vec3d(inter)};
        }, py::arg("s_s"), py::arg("s_e"), py::arg("plane_p"), py::arg("plane_n"),
        "3D线段与平面求交。返回 (是否相交, 交点[x,y,z])。");

    m.def("HGP_3D_Projection_Point_Plane_C1",
        [](const std::vector<double>& p,
           const std::vector<double>& pp, const std::vector<double>& pn)
           -> std::vector<double> {
            return from_vec3d(HGP_3D_Projection_Point_Plane_C1(
                to_vec3d(p), to_vec3d(pp), to_vec3d(pn)));
        }, py::arg("p"), py::arg("plane_p"), py::arg("plane_n"),
        "3D点投影到平面（C1，点+法向定义平面）。返回投影点 [x,y,z]。");

    m.def("HGP_3D_One_Triangle_Area",
        [](const std::vector<double>& v0,
           const std::vector<double>& v1,
           const std::vector<double>& v2) -> double {
            return HGP_3D_One_Triangle_Area(to_vec3d(v0), to_vec3d(v1), to_vec3d(v2));
        }, py::arg("v_0"), py::arg("v_1"), py::arg("v_2"), "单个3D三角形面积。");

    m.def("HGP_3D_Mesh_Curvature_C1",
        [](const std::vector<std::vector<double>>& vecs,
           const std::vector<int>& f0,
           const std::vector<int>& f1,
           const std::vector<int>& f2)
           -> std::pair<std::vector<double>, std::vector<double>> {
            Vector1d1 max_curs, min_curs;
            HGP_3D_Mesh_Curvature_C1(to_vec3d1(vecs), f0, f1, f2, max_curs, min_curs);
            return {max_curs, min_curs};
        }, py::arg("vecs"),
           py::arg("face_id_0"), py::arg("face_id_1"), py::arg("face_id_2"),
        "网格曲率C1（分量传参）。返回 (max_curs, min_curs)。");
        // ════════════════════════════════════════════════════
    // ── CSG: Primitive Generators
    // ════════════════════════════════════════════════════

    // HGP_Mesh_Make_Box
    // 返回 (vecs, fi0, fi1, fi2)
    m.def("HGP_Mesh_Make_Box",
        [](double cx, double cy, double cz,
           double wx, double wy, double wz)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>,
                         std::vector<int>,
                         std::vector<int>> {
            Vector3d1 vecs;
            Vector1i1 fi0, fi1, fi2;
            HGP_Mesh_Make_Box(cx, cy, cz, wx, wy, wz, vecs, fi0, fi1, fi2);
            return {from_vec3d1(vecs), fi0, fi1, fi2};
        },
        py::arg("cx"), py::arg("cy"), py::arg("cz"),
        py::arg("wx"), py::arg("wy"), py::arg("wz"),
        "生成长方体网格。(cx,cy,cz)=中心，(wx,wy,wz)=各轴半长。\n"
        "返回 (vecs, fi0, fi1, fi2)。");

    // HGP_Mesh_Make_Cylinder
    m.def("HGP_Mesh_Make_Cylinder",
        [](double cx, double cy, double cz,
           double radius, double height, int segments)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>,
                         std::vector<int>,
                         std::vector<int>> {
            Vector3d1 vecs;
            Vector1i1 fi0, fi1, fi2;
            HGP_Mesh_Make_Cylinder(cx, cy, cz, radius, height, segments,
                                   vecs, fi0, fi1, fi2);
            return {from_vec3d1(vecs), fi0, fi1, fi2};
        },
        py::arg("cx"), py::arg("cy"), py::arg("cz"),
        py::arg("radius"), py::arg("height"),
        py::arg("segments") = 32,
        "生成圆柱体网格。(cx,cy,cz)=中轴中点，segments 建议 >=32。\n"
        "返回 (vecs, fi0, fi1, fi2)。");

    // HGP_Mesh_Make_Sphere
    m.def("HGP_Mesh_Make_Sphere",
        [](double cx, double cy, double cz,
           double radius, int rings, int segments)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>,
                         std::vector<int>,
                         std::vector<int>> {
            Vector3d1 vecs;
            Vector1i1 fi0, fi1, fi2;
            HGP_Mesh_Make_Sphere(cx, cy, cz, radius, rings, segments,
                                 vecs, fi0, fi1, fi2);
            return {from_vec3d1(vecs), fi0, fi1, fi2};
        },
        py::arg("cx"), py::arg("cy"), py::arg("cz"),
        py::arg("radius"),
        py::arg("rings")    = 16,
        py::arg("segments") = 32,
        "生成 UV 球体网格。rings>=2，segments>=3。\n"
        "返回 (vecs, fi0, fi1, fi2)。");

    // HGP_Mesh_Make_Cone
    m.def("HGP_Mesh_Make_Cone",
        [](double cx, double cy, double cz,
           double radius, double height, int segments)
           -> std::tuple<std::vector<std::vector<double>>,
                         std::vector<int>,
                         std::vector<int>,
                         std::vector<int>> {
            Vector3d1 vecs;
            Vector1i1 fi0, fi1, fi2;
            HGP_Mesh_Make_Cone(cx, cy, cz, radius, height, segments,
                               vecs, fi0, fi1, fi2);
            return {from_vec3d1(vecs), fi0, fi1, fi2};
        },
        py::arg("cx"), py::arg("cy"), py::arg("cz"),
        py::arg("radius"), py::arg("height"),
        py::arg("segments") = 32,
        "生成圆锥体网格。顶点在 (cx, cy, cz+height)，底面圆心在 (cx,cy,cz)。\n"
        "返回 (vecs, fi0, fi1, fi2)。");

    // ════════════════════════════════════════════════════
    // ── CSG: Boolean Operations
    // ════════════════════════════════════════════════════

    // HGP_Mesh_CSG_Union  —  A ∪ B
    m.def("HGP_Mesh_CSG_Union",
        [](const std::vector<std::vector<double>>& vecs_a,
           const std::vector<int>& fi0_a,
           const std::vector<int>& fi1_a,
           const std::vector<int>& fi2_a,
           const std::vector<std::vector<double>>& vecs_b,
           const std::vector<int>& fi0_b,
           const std::vector<int>& fi1_b,
           const std::vector<int>& fi2_b)
           -> std::tuple<bool,
                         std::vector<std::vector<double>>,
                         std::vector<int>,
                         std::vector<int>,
                         std::vector<int>> {
            Vector3d1 vecs_out;
            Vector1i1 fi0_out, fi1_out, fi2_out;
            bool ok = HGP_Mesh_CSG_Union(
                to_vec3d1(vecs_a), fi0_a, fi1_a, fi2_a,
                to_vec3d1(vecs_b), fi0_b, fi1_b, fi2_b,
                vecs_out, fi0_out, fi1_out, fi2_out);
            return {ok, from_vec3d1(vecs_out), fi0_out, fi1_out, fi2_out};
        },
        py::arg("vecs_a"), py::arg("fi0_a"), py::arg("fi1_a"), py::arg("fi2_a"),
        py::arg("vecs_b"), py::arg("fi0_b"), py::arg("fi1_b"), py::arg("fi2_b"),
        "CSG 并集 A ∪ B。输入网格须为封闭流形。\n"
        "返回 (ok, vecs, fi0, fi1, fi2)。ok=False 时结果无效。");

    // HGP_Mesh_CSG_Difference  —  A - B
    m.def("HGP_Mesh_CSG_Difference",
        [](const std::vector<std::vector<double>>& vecs_a,
           const std::vector<int>& fi0_a,
           const std::vector<int>& fi1_a,
           const std::vector<int>& fi2_a,
           const std::vector<std::vector<double>>& vecs_b,
           const std::vector<int>& fi0_b,
           const std::vector<int>& fi1_b,
           const std::vector<int>& fi2_b)
           -> std::tuple<bool,
                         std::vector<std::vector<double>>,
                         std::vector<int>,
                         std::vector<int>,
                         std::vector<int>> {
            Vector3d1 vecs_out;
            Vector1i1 fi0_out, fi1_out, fi2_out;
            bool ok = HGP_Mesh_CSG_Difference(
                to_vec3d1(vecs_a), fi0_a, fi1_a, fi2_a,
                to_vec3d1(vecs_b), fi0_b, fi1_b, fi2_b,
                vecs_out, fi0_out, fi1_out, fi2_out);
            return {ok, from_vec3d1(vecs_out), fi0_out, fi1_out, fi2_out};
        },
        py::arg("vecs_a"), py::arg("fi0_a"), py::arg("fi1_a"), py::arg("fi2_a"),
        py::arg("vecs_b"), py::arg("fi0_b"), py::arg("fi1_b"), py::arg("fi2_b"),
        "CSG 差集 A - B（从 A 中挖去 B）。输入网格须为封闭流形。\n"
        "返回 (ok, vecs, fi0, fi1, fi2)。");

    // HGP_Mesh_CSG_Intersection  —  A ∩ B
    m.def("HGP_Mesh_CSG_Intersection",
        [](const std::vector<std::vector<double>>& vecs_a,
           const std::vector<int>& fi0_a,
           const std::vector<int>& fi1_a,
           const std::vector<int>& fi2_a,
           const std::vector<std::vector<double>>& vecs_b,
           const std::vector<int>& fi0_b,
           const std::vector<int>& fi1_b,
           const std::vector<int>& fi2_b)
           -> std::tuple<bool,
                         std::vector<std::vector<double>>,
                         std::vector<int>,
                         std::vector<int>,
                         std::vector<int>> {
            Vector3d1 vecs_out;
            Vector1i1 fi0_out, fi1_out, fi2_out;
            bool ok = HGP_Mesh_CSG_Intersection(
                to_vec3d1(vecs_a), fi0_a, fi1_a, fi2_a,
                to_vec3d1(vecs_b), fi0_b, fi1_b, fi2_b,
                vecs_out, fi0_out, fi1_out, fi2_out);
            return {ok, from_vec3d1(vecs_out), fi0_out, fi1_out, fi2_out};
        },
        py::arg("vecs_a"), py::arg("fi0_a"), py::arg("fi1_a"), py::arg("fi2_a"),
        py::arg("vecs_b"), py::arg("fi0_b"), py::arg("fi1_b"), py::arg("fi2_b"),
        "CSG 交集 A ∩ B。输入网格须为封闭流形。\n"
        "返回 (ok, vecs, fi0, fi1, fi2)。");    
}