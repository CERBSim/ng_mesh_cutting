#include <pybind11/pybind11.h>

#include <meshing.hpp>

using namespace netgen;

void CutMesh(Mesh& mesh, FlatArray<double> values, int domain_to_split,
             int dom_pos, int dom_neg, int bc_number) {
  mesh.UpdateTopology();
  const auto& topology = mesh.GetTopology();

  constexpr int8_t P_NEG = -1, P_POS = 1, P_ZERO = 0;
  constexpr int8_t P_POS_NEG = P_NEG * P_POS;

  Array<int8_t, PointIndex> signs(mesh.GetNP());
  constexpr double eps = 1e-8;
  for (PointIndex pi : Range(mesh.Points())) {
    if (fabs(values[pi - 1]) < eps)
      signs[pi] = P_ZERO;
    else if (values[pi - 1] < 0)
      signs[pi] = P_NEG;
    else
      signs[pi] = P_POS;
    // cout << pi << '\t' << ((int)(signs[pi])) << endl;
  }

  auto fd = FaceDescriptor(dom_neg - 1, 1, 0, 0);
  fd.SetBCProperty(dom_neg);
  mesh.AddFaceDescriptor(fd);
  auto elementsonnode = mesh.CreatePoint2SurfaceElementTable();
  Array<std::tuple<PointIndex, PointIndex>> edges;
  BuildEdgeList(mesh, elementsonnode, edges);

  std::map<std::tuple<PointIndex, PointIndex>, PointIndex> point_map;

  for (auto [pi1, pi2] : edges) {
    if (pi2 < pi1) swap(pi1, pi2);
    // cout << "check " << pi1 << ' ' << pi2 << endl;
    if (signs[pi1] * signs[pi2] == P_POS_NEG) {
      auto val1 = values[pi1 - 1];
      auto val2 = values[pi2 - 1];
      // cout << "\tadd point" << pi1 << ' ' << pi2 << endl;
      auto& p1 = mesh.Point(pi1);
      auto& p2 = mesh.Point(pi2);
      auto t = val1 / (val1 - val2);
      auto p = p1 + t * (p2 - p1);
      auto pi_new = mesh.AddPoint(p);
      point_map[{pi1, pi2}] = pi_new;
      // cout << "insert " << pi1 << ", " << pi2 << endl;
      auto segs = topology.GetVertexSegments(pi1);
      SegmentIndex segi = -1;
      for (auto si : segs) {
        auto& s = mesh.LineSegment(si);
        if (s[0] == pi2 || s[1] == pi2) {
          segi = si;
          break;
        }
      }
      if (segi != -1) {
        auto seg_old = mesh.LineSegment(segi);
        for (auto i : Range(2)) {
          Segment seg_new = seg_old;
          seg_new[i] = pi_new;
          // if(signs[seg_new[1-i]] > 0) {
          //   seg_new.edgenr = 5;
          //   seg_new.si = 5;
          // }
          // cout << "split segment " << seg_new << endl;
          mesh.AddSegment(seg_new);
        }
        mesh.DeleteSegment(segi + 1);
      }
    }
  }

  std::map<PointIndex, ArrayMem<SegmentIndex, 2>> point_to_segments;

  auto nsel = mesh.SurfaceElements().Size();
  for (SurfaceElementIndex sei : Range(nsel)) {
    auto sel = mesh[sei];
    ArrayMem<PointIndex, 4> p_pos, p_neg, p_zero;
    for (auto pi : sel.PNums()) {
      if (signs[pi] == P_NEG)
        p_neg.Append(pi);
      else if (signs[pi] == P_POS)
        p_pos.Append(pi);
      else
        p_zero.Append(pi);
    }
    if (p_pos.Size() > 0 && p_neg.Size() > 0) {
      if (sel.GetIndex() != domain_to_split)
        throw Exception("domain_to_split is not correct");
      if (p_zero.Size() == 0) {
        auto p0 = p_pos.Size() == 1 ? p_pos[0] : p_neg[0];
        auto p1 = p_pos.Size() == 1 ? p_neg[0] : p_pos[0];
        auto p2 = p_pos.Size() == 1 ? p_neg[1] : p_pos[1];

        bool need_invert = sel.PNums().Pos(p0) == 1;

        // cout << "get " << min(p0, p1) << ", " << max(p0, p1) << endl;
        auto p01 = point_map.at({min(p0, p1), max(p0, p1)});
        // cout << "get " << min(p0, p2) << ", " << max(p0, p2) << endl;
        auto p02 = point_map.at({min(p0, p2), max(p0, p2)});
        Element2d el;
        el = sel;
        el.PNums() = ArrayMem<PointIndex, 3>{p0, p01, p02};
        int dom = signs[p0] == P_NEG ? dom_neg : dom_pos;
        el.SetIndex(dom);
        if (need_invert) el.Invert();
        mesh.AddSurfaceElement(el);

        el.PNums() = ArrayMem<PointIndex, 3>{p01, p1, p2};
        dom = dom_neg + dom_pos - dom;
        el.SetIndex(dom);
        if (need_invert) el.Invert();
        mesh.AddSurfaceElement(el);

        el.PNums() = ArrayMem<PointIndex, 3>{p2, p02, p01};
        el.SetIndex(dom);
        if (need_invert) el.Invert();
        mesh.AddSurfaceElement(el);

        Segment seg;
        seg[0] = p01;
        seg[1] = p02;
        seg.edgenr = bc_number;
        seg.si = bc_number;
        // SetBCName(5, "split");

        // cout << "points " << p0 << ", " << p1 << ", " << p2 << ", " << p01 <<
        // ", " << p02 << endl; cout << "seg2 " << p0 << ", " << p1 << ", " <<
        // p2 << ", " << p01 << ", " << p02 << endl;
        auto segi = mesh.AddSegment(seg);
        point_to_segments[p01].Append(segi);
        point_to_segments[p02].Append(segi);
      }
      if (p_zero.Size() == 1) {
        auto p0 = p_zero[0];
        auto p1 = p_neg[0];
        auto p2 = p_pos[0];

        bool need_invert = sel.PNums().Pos(p0) == 1;

        auto p12 = point_map.at({min(p1, p2), max(p1, p2)});

        Element2d el;
        el = sel;
        el.PNums() = ArrayMem<PointIndex, 3>{p0, p1, p12};
        el.SetIndex(dom_neg);
        if (need_invert) el.Invert();
        mesh.AddSurfaceElement(el);

        el.PNums() = ArrayMem<PointIndex, 3>{p0, p12, p2};
        el.SetIndex(dom_pos);
        if (need_invert) el.Invert();
        mesh.AddSurfaceElement(el);

        Segment seg;
        seg[0] = p0;
        seg[1] = p12;
        seg.edgenr = bc_number;
        seg.si = bc_number;
        // SetBCName(5, "split");
        // cout << "points " << p0 << ", " << p1 << ", " << p2 << ", " << p12 <<
        // endl; cout << "add segment " << seg << endl;
        auto segi = mesh.AddSegment(seg);
        point_to_segments[p0].Append(segi);
        point_to_segments[p12].Append(segi);
      }
      if (p_zero.Size() == 2) {
        cout << "impossible!" << endl;
        throw Exception();
      }

      mesh[sei].Delete();
    } else {
      auto& sel = mesh[sei];
      if (sel.GetIndex() == domain_to_split)
        for (auto pi : sel.PNums()) {
          if (signs[pi] == P_POS) sel.SetIndex(dom_pos);
          if (signs[pi] == P_NEG) sel.SetIndex(dom_neg);
        }
    }
  }
  // mesh.Save("out.vol");

  mesh.Compress();
  mesh.SetNextMajorTimeStamp();
  mesh.RebuildSurfaceElementLists();
  mesh.CalcSurfacesOfNode();
  mesh.UpdateTopology();

  // mesh.Save("mesh.vol");
}

PYBIND11_MODULE(ng_mesh_cutting, m) { m.def("CutMesh", &CutMesh); }
