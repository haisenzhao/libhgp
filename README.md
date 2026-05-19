<div align="center">

# LibHGP

**A Lightweight Geometry Engine for Digital Manufacturing**

A foundational algorithm library for integrating and industrializing laboratory geometric algorithm assets, supporting **C++ / Python / Web** three calling interfaces.

</div>

## Table of Contents

- [Introduction](#introduction)
- [Architecture](#architecture)
- [Recommended Environment](#recommended-environment)
- [Quick Start](#quick-start)
- [License](#license)

## Introduction

LibHGP is a geometry processing and computational geometry algorithm library that serves as the foundation for building an interactive visualization platform with good maintainability, scalability, and reusability. It supports three usage modes:

- **C++ Native Calls**: Suitable for high-performance algorithm development and system integration
- **Python Calls**: Suitable for research experiments, prototyping, and teaching
- **Web Calls**: Suitable for frontend display, interactive visualization, and online demonstrations

## Architecture

The project can be divided into the following layers:

### 1. Core Algorithm Layer
Responsible for implementing geometry processing, computational geometry, and related fundamental algorithm modules, including:

- 2D geometry algorithms
- 3D geometry algorithms
- Mesh processing
- Convex hull and spatial indexing
- Other laboratory-accumulated algorithm modules

### 2. Interface Wrapper Layer
Responsible for providing unified wrappers for underlying algorithms, exposing stable and consistent calling interfaces, including:

- C++ native interface
- Python binding interface (e.g., pybind11)
- Web calling interface (can support FastAPI / WebAssembly in the future)

### 3. Application and Presentation Layer
Targeting specific application scenarios, including:

- Research prototyping
- Teaching demonstrations
- Algorithm visualization
- Web-based interactive display

## Recommended Environment
- C++17 
- Visual Studio 2019

## Quick Start

### C++ Users

1. Copy all files from `libhgp/C++/` into your project directory
2. Include the header file: `#include "libhgp.h"`
3. Declare the namespace: `using namespace libhgp;`
4. Use the following code for a quick test:
```cpp
#include "iostream"
#include"libhgp.h"
using namespace std;
using namespace libhgp;
int main(int argc, char* argv[])
{
	char dirBuf[260] = { 0 };
	char userBuf[260] = { 0 };
	auto a = PL().HGP_2D_Distance_Point_Point_C(Vector2d(0, 0), Vector2d(1, 1));
	PL().LGP_WinGetCurDirectory_C(dirBuf, sizeof(dirBuf));
	PL().LGP_WinGetUserName_C(userBuf, sizeof(userBuf));
	std::cerr << a << std::endl;
	cerr << "WinGetCurDirectory: " << dirBuf << endl;
	cerr << "WinGetUserName: " << userBuf << endl;
	system("pause");
	return 0;
}
```
### Python Users

1. Copy all files from `libhgp/Python/` into your project directory
2. Import the module in Python: `import hgp_py`
3. Use the following code for a quick test:
```python
import hgp_py

# ── Test 1: 2D point-to-point distance ──────────────────────────────
p1, p2 = [0.0, 0.0], [3.0, 4.0]
dist = hgp_py.HGP_2D_Distance_Point_Point(p1, p2)
assert abs(dist - 5.0) < 1e-9
print(f"[PASS] Test 1: 2D point distance (0,0) (3,4) \n HGP_2D_Distance_Point_Point({p1}, {p2}) = {dist} (expected 5.0)")

# ── Test 2: Point-in-polygon test ────────────────────────────
polygon = [[0.0,0.0],[4.0,0.0],[4.0,4.0],[0.0,4.0]]
inside_pt, outside_pt = [2.0, 2.0], [5.0, 5.0]
assert hgp_py.HGP_2D_Location_Point_Polygon(inside_pt, polygon) == True
assert hgp_py.HGP_2D_Location_Point_Polygon(outside_pt, polygon) == False
print(f"[PASS] Test 2: Point-in-polygon test Polygon[(0,0),(4,0),(4,4),(0,4)] \n HGP_2D_Location_Point_Polygon({inside_pt}) = True")
print(f" HGP_2D_Location_Point_Polygon({outside_pt}) = False")

# ── Test 3: Polygon regular sampling ──────────────────────────────
pts = hgp_py.HGP_2D_Polygon_Regular_Sampling_C1(polygon, 0.1)
print(f"[PASS] Test 3: Polygon regular sampling Polygon[(0,0),(4,0),(4,4),(0,4)] spacing=0.1 \n HGP_2D_Polygon_Regular_Sampling_C1() → generated {len(pts)} sample points")

# ── Test 4: Line segment intersection ────────────────────────────────────
seg1_p1, seg1_p2 = [0.0,0.0],[2.0,2.0]
seg2_p1, seg2_p2 = [0.0,2.0],[2.0,0.0]
hit, inter = hgp_py.HGP_2D_Intersection_Segment_Segment(
    seg1_p1, seg1_p2, seg2_p1, seg2_p2)
assert hit == True
assert abs(inter[0] - 1.0) < 1e-9 and abs(inter[1] - 1.0) < 1e-9
print(f"[PASS] Test 4: Line segment intersection Segment A[(0,0)→(2,2)] Segment B[(0,2)→(2,0)] \n HGP_2D_Intersection_Segment_Segment() = {inter} (expected [1.0,1.0])")

# ── Test 5: Polygon area ──────────────────────────────────
area = hgp_py.HGP_2D_Polygon_Area(polygon)
print(f"[PASS] Test 5: Polygon area Polygon[(0,0),(4,0),(4,4),(0,4)] \n HGP_2D_Polygon_Area() = {area} (expected 16.0)")

# ── Test 6: 3D point-to-point distance ────────────────────────────────
p3, p4 = [0.0,0.0,0.0], [1.0,1.0,1.0]
d3 = hgp_py.HGP_3D_Distance_Point_Point(p3, p4)
expected = 3**0.5
print(f"[PASS] Test 6: 3D point distance (0,0,0) (1,1,1) \n HGP_3D_Distance_Point_Point({p3}, {p4}) = {d3:.6f} (expected √3 ≈ {expected:.6f})")

print("\n✅ All tests passed.")
"
```
### Web Users

1. Copy all files from `libhgp/Web/` into your project directory
2. Navigate to the server build directory: `cd backend`
3. Start the local server: `python -m uvicorn app:app --reload --port 8000`
4. Access the local server at: `http://127.0.0.1:8000/`

## License

All rights about the program are reserved by the authors of this project. The programs can only be used for research purpose. In no event shall the author be liable to any party for direct, indirect, special, inc
