#include<phat/boundary_matrix.h>
#include<phat/compute_persistence_pairs.h>

#include "persistence_module.h"

namespace md {
    PersistenceModule::PersistenceModule(const std::string& /*fname*/) // read from file
            :
            generators_(),
            relations_()
    {
    }

    Diagram PersistenceModule::slice_diagram(const DualPoint& /*line*/)
    {
        //Vector2D b_of_line(L.b, -L.b);
        //for (int i = 0; i<(int) F.size(); i++) {
        //    Simplex_in_2D_filtration& curr_simplex = F[i];
        //    Vector2D proj = push(curr_simplex.pos, L);

        //    curr_simplex.v = L_2(proj - b_of_line);
        //    //std::cout << proj << std::endl;
        //    //std::cout << "v=" << curr_simplex.v << std::endl;
        //}
        //std::sort(F.begin(), F.end(), sort_functor);
        //std::map<Index, Index> index_map;
        //for (Index i = 0; i<(int) F.size(); i++) {
        //    index_map[F[i].index] = i;
        //    //std::cout << F[i].index << " -> " << i << std::endl;
        //}
        //phat::boundary_matrix<> phat_matrix;
        //phat_matrix.set_num_cols(F.size());
        //std::vector<Index> bd_in_slice_filtration;
        //for (Index i = 0; i<(int) F.size(); i++) {
        //    phat_matrix.set_dim(i, F[i].dim);
        //    bd_in_slice_filtration.clear();
        //    //std::cout << "new col" << i << std::endl;
        //    for (int j = 0; j<(int) F[i].bd.size(); j++) {
        //        // F[i] contains the indices of its facet wrt to the
        //        // original filtration. We have to express it, however,
        //        // wrt to the filtration along the slice. That is why
        //        // we need the index_map
        //        //std::cout << "Found " << F[i].bd[j] << ", returning " << index_map[F[i].bd[j]] << std::endl;
        //        bd_in_slice_filtration.push_back(index_map[F[i].bd[j]]);
        //    }
        //    std::sort(bd_in_slice_filtration.begin(), bd_in_slice_filtration.end());
        //    [>
        //    for(int j=0;j<bd_in_slice_filtration.size();j++) {
        //      std::cout << bd_in_slice_filtration[j] << " ";
        //    }
        //    std::cout << std::endl;
        //    */
        //    phat_matrix.set_col(i, bd_in_slice_filtration);

        //}
        //phat::persistence_pairs phat_persistence_pairs;
        //phat::compute_persistence_pairs<phat::twist_reduction>(phat_persistence_pairs, phat_matrix);

        //dgm.clear();
        //for (long i = 0; i<(long) phat_persistence_pairs.get_num_pairs(); i++) {
        //    std::pair<phat::index, phat::index> new_pair = phat_persistence_pairs.get_pair(i);
        //    double birth = F[new_pair.first].v;
        //    double death = F[new_pair.second].v;
        //    if (birth!=death) {
        //        dgm.push_back(std::make_pair(birth, death));
        //    }
        //}

        //[>
        //std::cout << "Done, created diagram: " << std::endl;
        //for(int i=0;i<(int)dgm.size();i++) {
        //  std::cout << dgm[i].first << " " << dgm[i].second << std::endl;
        //}
        //*/
        return Diagram();

    }

    PersistenceModule::Box PersistenceModule::bounding_box() const
    {
        Real ll_x = std::numeric_limits<Real>::max();
        Real ll_y = std::numeric_limits<Real>::max();
        Real ur_x = -std::numeric_limits<Real>::max();
        Real ur_y = -std::numeric_limits<Real>::max();

        for(const auto& gen : generators_) {
            ll_x = std::min(gen.x, ll_x);
            ll_y = std::min(gen.y, ll_y);
            ur_x = std::max(gen.x, ur_x);
            ur_y = std::max(gen.y, ur_y);
        }

        for(const auto& rel : relations_) {

            ll_x = std::min(rel.get_x(), ll_x);
            ll_y = std::min(rel.get_y(), ll_y);

            ur_x = std::max(rel.get_x(), ur_x);
            ur_y = std::max(rel.get_y(), ur_y);
        }

        return Box(Point(ll_x, ll_y), Point(ur_x, ur_y));
    }
}
