#ifndef TREE_NODE_H
#define TREE_NODE_H

#include <cfloat>
#include <stdint.h>

namespace RoutingKit {
    class tree_node {
    private:
        unsigned id_;
        unsigned parent_id_;
        unsigned g_;
        bool status_;  // open or closed

    public:
        // tree_node() = default;
        tree_node() { g_ = INT32_MAX; }

        tree_node(unsigned id, unsigned parent_id, unsigned g)
                : id_{id}, parent_id_{parent_id}, g_{g} {
            status_ = 0;
        }

        ~tree_node() {}

        // Compile error, why? Forward declaration
        // tree_node(const tree_node&) = delete;

        // tree_node& operator=(const tree_node&) = delete;

        inline uint32_t get_id() const { return id_; }

        inline void set_id(uint32_t id) { id_ = id; }

        inline uint32_t get_parent() const { return parent_id_; }

        inline void set_parent(uint32_t parent_id) { parent_id_ = parent_id; }

        inline double get_g() const { return g_; }

        inline void set_g(double g) { g_ = g; }

        inline bool get_expanded() const { return status_; }

        inline void set_expanded(bool expanded)  { status_ = expanded; }
    };
}  // namespace warthog
#endif
