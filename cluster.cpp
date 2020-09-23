#include <iostream>
#include <cmath>
#include <random>
#include <limits>
#include <cmath>
#include <utility>
#include <algorithm>
#include <set>

#define SQRT2PI (2.506628274631001)
#define SQRTPI (1.772453851)
#define SQRT2 (1.414213562373095)

using std::max;

class cluster;

std::ostream & operator<<(std::ostream &, cluster const &);

struct dist {
    double mean;
    double m2;
    int count;

    inline double population_var() const {
        return m2 / count;
    }

    inline double sample_var() const {
        // have to apply bessels correction
        if (count == 1) {
            // return 0;
            return std::numeric_limits<double>::max();
        }
        return m2 / (count - 1);
    }

    inline double sample_std() const {
        return std::sqrt(sample_var());
    }
    inline double population_std() const {
        return std::sqrt(population_var());
    }

    void update(double sample) {
        count++;
        double delta = sample - mean;
        mean += delta / count;
        double delta2 = sample - mean;
        m2 += delta * delta2;
    }

    double sample_density(double sample) const {
        // if we have a single sample, the density everywhere is 1/sqrt(2 * PI) because it has infinite variance
        if (count == 1) {
            return 1 / SQRT2PI;
        }

        double delta = sample - mean;
        // should we include the count or not?
        // it's going to favor larger clusters, so maybe we should ignore it
        return std::exp(-delta * delta / sample_var() / 2.) / SQRT2PI / sample_std();
    }

    double population_density(double sample) const {
        double delta = sample - mean;
        // should we include the count or not?
        // it's going to favor larger clusters, so maybe we should ignore it
        return std::exp(-delta * delta / population_var() / 2.) / SQRT2PI / sample_std();
    }



    dist() : mean(0), m2(0), count(0) { }
    dist(double mean) : mean(mean), m2(0), count(1) { }
    dist(double mean, double stddev, int count) : mean(mean), m2(stddev * stddev * count), count(count) { }
    dist(dist const & d) : mean(d.mean), m2(d.m2), count(d.count) { }
    dist & operator=(dist const & d) {
        if(this != &d) {
            mean = d.mean;
            m2 = d.m2;
            count = d.count;
        }
        return *this;
    }
    dist(dist const & A, dist const & B) {
        count = A.count + B.count;
        double a = (double)A.count / (double)count;
        double b = 1. - a;
        double delta = A.mean - B.mean;

        mean = a * A.mean + b * B.mean;
        double var = a * A.population_var() + b * B.population_var() + a * b * delta * delta;
        m2 = var * count;
    }

};

std::ostream& operator<<(std::ostream& os, dist const & d) {
    return os << "[" << d.mean << ";" << d.population_std() << "#" << d.count << "]";
}

dist mix(dist const & A, dist const & B) {
    double c = A.count + B.count;
    double a = (double)A.count / c;
    double b = 1. - a;

    double delta = A.mean - B.mean;

    double mean = a * A.mean + b * B.mean;
    double var  = a * A.population_var() + b * B.population_var() + a * b * delta * delta;

    dist r;
    r.mean = mean;
    r.m2 = var * c;
    r.count = (int)c;
    
    // std::cout << "mix: " << A << " " << B << " -> " << r << std::endl;
    // std::cout << "a: " << a << " b: " << b << " delta: " << delta << std::endl;

    return r;
}

// returns true if the mixture of A and B would create a distribution with a variance greater than both A and B
bool can_mix(dist const & A, dist const & B) {
    // if (A.count == 1 || B.count == 1) return true;
    // let's do this the dumb way first:

    dist m = mix(A, B);
    return m.population_var() > A.population_var() && m.population_var() > B.population_var();

    // double c = A.count + B.count;
    // double Avar = A.var();
    // double Bvar = B.var();
    // double alpha, deltaVar;
    // if (Avar > Bvar) {
    //     alpha = (double)A.count / c;
    //     deltaVar = Avar - Bvar;
    // } else {
    //     alpha = (double)B.count / c;
    //     deltaVar = Bvar - Avar;
    // }
    // return alpha * (A.mean - B.mean) * (A.mean - B.mean) > deltaVar;
}

double gaussian(double mean, double variance, double x) {
    double d = x - mean;
    return std::exp(-d * d / 2 / variance) / SQRT2PI / std::sqrt(variance);
}

/*
returns a number between 0 if the distributions are identical, 1 if they do not overlap at all
*/
double comp_dist(dist const & left, dist const & right) {
    if (left.population_var() == 0 || right.population_var() == 0) {
        return 1;
    } 

    double d = left.mean - right.mean;
    double v = left.population_var() + right.population_var();
    double s = left.population_std() + right.population_std();
    double sa = left.population_std();
    double sb = right.population_std();

    // this should never be less than zero, if it is, it must be zero
    double x = 1 - 2. * SQRT2 * sa * sb * std::exp(- d * d / 2. / v) / std::sqrt(v) / s;

    return (x < 0) ? 0 : (x > 1) ? 1 : x;
}

class cluster {
    friend std::ostream & operator<<(std::ostream &, cluster const &);
protected:
    struct node {
        node * left, * right;
        dist d;

        node(double m) : 
            left(nullptr),
            right(nullptr),
            d(m)
        { }

        node(node * left, node * right) :
            left(left), right(right),
            d(mix(left->d, right->d))
        { 
            if (left->d.population_var() > right->d.population_var()) {
                std::swap(left, right);
            }
        }

        node(node * left, node * right, dist d) 
          : left(left), right(right), d(d)
        { 
            if (left->d.population_var() > right->d.population_var()) {
                std::swap(left, right);
            }
        }

        double error2() const {
            if (left == nullptr) {
                return 0;
            }

            double a = (double)left->d.count / (double)d.count;
            double b = 1. - a;

            // double va = left->d.sample_var();
            // double vb = right->d.sample_var();
            // double vc = d.sample_var();
            double sa = left->d.sample_std();
            double sb = right->d.sample_std();
            double sc = d.sample_std();
            double va = left->d.sample_var();
            double vb = right->d.sample_var();
            double vc = d.sample_var();
            double ma = left->d.mean;
            double mb = right->d.mean;
            double mc = d.mean;

            if (sa == 0 || sb == 0 || sc == 0) {
                return std::numeric_limits<double>::max();
                //throw std::logic_error("stddev zero");
            }

            double ret = 0.;

            ret += a * a / sa;
            ret += b * b / sb;
            ret += 1 / sc;
            ret = ret / 2 / SQRTPI;

            ret += 2 * a * b * gaussian(ma, va + vb, mb);
            ret -= 2 * a * gaussian(ma, va + vc, mc);
            ret -= 2 * b * gaussian(mb, vb + vc, mc);

            return ret;
        }
    };

private:
    node * root;
    std::minstd_rand gen;
    std::uniform_real_distribution<double> uniform;

    void insert_helper3(node ** from, node ** sample) {
        typedef void (cluster::*RecurseType)(node**,node**);
        RecurseType recurse = cluster::insert_helper3;

        if (*sample == nullptr) return;

        if (*from == nullptr) {
            *from = *sample;
            return;
        }

        node * f = *from;
        node * s = *sample;

        if (f->left == nullptr) {
            dist m = mix(f->d, s->d);

            if (m.population_var() < max(f->d.population_var(), s->d.population_var())) {
                // ok, first problem. We are at a leaf in the spot where
                // we feel it's best to insert, but we can't mix with the
                // sample.  What should we do?

                // option 1: bubble up the one with the higher variance
                //   honestly I can't think of another option

                if (f->d.population_var() > s->d.population_var()) {
                    *from = s;
                    *sample = f;
                }
            } else {
                *from = new node(f, s);
                *sample = nullptr;
            }
            return;
        }

        // we are not a leaf, so we could insert here or keep going
        // how do we decide?
        // I don't think we should ever insert here until we've tried to insert
        // into one of the child nodes.

        auto l = comp_dist(f->left->d, s->d);
        auto r = comp_dist(f->right->d, s->d);
        bool insertLeft = true;

        // less is more similar
        if (l < r) {
            insert_helper3(&f->left, sample);
        } else {
            insert_helper3(&f->right, sample);
            insertLeft = false;
        }

        // keep the lower variance to the left
        if (f->left->d.population_var() > f->right->d.population_var()) {
            std::swap(f->left, f->right);
        }

        s = *sample;

        // first check to see if we can mix the two new children
        auto m = mix(f->left->d, f->right->d);

        // can left ad right mix?
        if (m.population_var() < max(f->left->d.population_var(), f->right->d.population_var())) {
            // ok, we can't mix left and right  let's check the sample

            if (s == nullptr) {
                // no bubble, let's just bubble up the larger population variance
                if (f->left->d.population_var() < f->right->d.population_var()) {
                    *from = f->left;
                    *sample = f->right;
                } else {
                    *from = f->right;
                    *sample = f->left;
                }
                return;
            }

            // here is the trickiest bit, we have a bubble, and we can't mix.
            // first lets see if we can mix bubble with either left or right:

            if (insertLeft) {
                // first we'll check to see if we can insert s here:
                auto n = mix(f->left->d, s->d);

                if (n.population_var() > max(f->left->d.population_var(), s->d.population_var())) {
                    // cool, just merge left with the bubble, and bubble right
                    *from = new node(f->left, s);
                    *sample = f->right;
                    delete f;
                    return;
                }

                // ok, that didn't work, so try right
                n = mix(f->right->d, s->d);

                if (n.population_var() > max(f->right->d.population_var(), s->d.population_var())) {
                    // splice in the bubble and bubble up the left one
                    *from = new node(s, f->right);
                    *sample = f->left;
                    delete f;
                    return;
                }

                // here's the really tricky part.  we have a bubble but can't mix the other one
                // so we have three distributions that don't mix in pairs.
                // s came from the left side, and the left side is not a leaf.

                // try rotating the tree (why do I think that will work?)
                // before we inserted mix(left,right) > max(left,right)
                // we also know that s.var > max(left->left.var, left->right.var)
                
            } else {
                dist n = mix(f->left->d, s->d);

                if (n.population_var() > max(f->left->d.population_var(), s->d.population_var())) {
                    // splice in the bubble and bubble up the right one;
                    *from = new node(s, f->left);
                    *sample = f->right;
                    delete f;
                    return;
                }
            }


        } else {
            f->d = m;
        }

        // do we have a bubble?
        if (s != nullptr) {
            // can we mix this?
            auto m = mix(f->d, s->d);

            if (m.population_var() < max(f->d.population_var(), s->d.population_var())) {
                // we can't mix, let it bubble
                return;
            } else {
                *from = new node(f, s);
                *sample = nullptr;
                // we are good here, return!
                return;
            }
        }



        // we can't mix, but left and right are good.  
        // can we rotate the tree to make it work?
    }

public:
    void insert3(double sample) {
        node * s = new node(sample);

        node * bub = nullptr;
        insert_helper3(&root, &s);

        // if (s != nullptr) {
        //     print_helper(std::cout, 0, root);
        //     print_helper(std::cout, 0, bub);

        //     throw std::logic_error("got a bubble at root!");
        // }
    }

private:

    void insert_helper2(node ** from, node * sample, node ** bubble) {
        node * f = *from;

        // is this a trailing edge?
        if (f == nullptr) {
            *from = sample;
            return;
        }

        node * b = nullptr;
        bool insertedLeft = true;

        // is this a leaf?
        if (f->left == nullptr) {
            *from = new node(f, sample);
            f = *from;
        } else {
            // which variance will increase more?
            dist ld = mix(sample->d, f->left->d);
            dist rd = mix(sample->d, f->right->d);

            if (ld.population_var() < rd.population_var()) {
                insert_helper2(&f->left, sample, &b);
            } else {
                insert_helper2(&f->right, sample, &b);
                insertedLeft = false;
            }
        }

        if (f->left != nullptr) {
            // now check our condition
            dist m = mix(f->left->d, f->right->d);

            if (m.population_var() < max(f->left->d.population_var(), f->right->d.population_var())) {
                if (b != nullptr) {
                    // can either one mix with b?
                    dist n = mix(b->d, f->left->d);
                    dist o = mix(b->d, f->right->d);

                    if(n.population_var() > max(b->d.population_var(), f->left->d.population_var())) {
                        // we can merge left with the bubble
                        *from = f->right;
                        b = new node(b, f->left);
                        delete f;
                        f = *from;
                    } else if (o.population_var() > max(b->d.population_var(), f->right->d.population_var())) {
                        // we can merge right with the bubble
                        *from = f->left;
                        b = new node(b, f->right);
                        delete f;
                        f = *from;
                    } else {
                        // std::cout << "can't fix and have a bubble: \n b:" << b->d << " l:" << f->left->d << " r:" << f->right->d << std::endl;
                        // std::cout << "lr:" << m << " bl:" << n << " br:" << o << std::endl;
                        // print_helper(std::cout, 0, b);
                        // print_helper(std::cout, 0, f);

                        node * nb = nullptr;

                        if (n.population_var() < o.population_var()) {
                            // try inserting f->left into b
                            *from = f->right;
                            insert_helper2(&b, f->left, &nb);

                            if (nb != nullptr) {
                                std::cout << "tried inserting left into bubble but a new bubble formed!" << std::endl;
                                throw std::logic_error("can't merge and have a bubble");
                            } 
                            delete f;
                            f = *from;
                        } else {
                            *from = f->left;
                            insert_helper2(&b, f->right, &nb);

                            if (nb != nullptr) {
                                std::cout << "tried inserting right into bubble but a new bubble formed!" << std::endl;
                                throw std::logic_error("can't merge and have a bubble");
                            }
                            delete f;
                            f = *from;
                        }

                    }
                } else {
                    if (f->left->d.population_var() < f->right->d.population_var()) {
                        b = f->right;
                        *from = f->left;
                        delete f;
                        f = *from;
                    } else {
                        b = f->left;
                        *from = f->right;
                        delete f;
                        f = *from;
                    }
                }
            } else {
                f->d = m;
            }
        }
        
        if (b != nullptr) {
            // can we insert the bubble here?
            dist m = mix(f->d, b->d);

            if(m.population_var() > max(f->d.population_var(), b->d.population_var())) {
                // we can insert it here!
                *from = new node(f, b);
                b = nullptr;
            } else {
                *bubble = b;
            }
        }

    }

    bool insert_helper(node ** from, node * sample, node ** bubble) {
        if (*from == nullptr) {
            *from = sample;
            return true;
        }

        // // can sample mix with *from?
        // if (!can_mix((*from)->d, sample->d)) {
        //     // if we can't mix, bubble up the one with the biggest variance
        //     std::cout << "can't mix nodes: " << (*from)->d << " " << sample->d << " -> " << mix((*from)->d, sample->d);
            
        //     if ((*from)->d.population_var() > sample->d.population_var()) {
        //         std::cout << " bubbling self.";
        //         *bubble = *from;
        //         *from = sample;
        //     } else {
        //         std::cout << " bubbling sample.";
        //         // why would this happen?
        //         *bubble = sample;
        //     }
        //     std::cout << std::endl;
        //     return;
        // }

        // is this a leaf
        if ((*from)->left == nullptr) {
            if (can_mix((*from)->d, sample->d)) {
                *from = new node(*from, sample);
                return true;
            } else {
                // if we can't mix, bubble up the one with the biggest variance
                std::cout << "can't mix nodes: " << (*from)->d << " " << sample->d << " -> " << mix((*from)->d, sample->d);
                
                if ((*from)->d.population_var() > sample->d.population_var()) {
                    std::cout << " bubbling self.";
                    *bubble = *from;
                    *from = sample;
                } else {
                    std::cout << " bubbling sample.";
                    // why would this happen?
                    *bubble = sample;
                }
                std::cout << std::endl;
                return true;
            }
        }

        // or else recurse
        // calculate the probability based on the density
        double leftVal = gaussian((*from)->left->d.mean, (*from)->left->d.population_var() + sample->d.population_var(), sample->d.mean);
        double rightVal = gaussian((*from)->right->d.mean, (*from)->right->d.population_var() + sample->d.population_var(), sample->d.mean);

        // flip a weighted coin
        double x = (leftVal + rightVal) * uniform(gen);

        node * bub = nullptr;
        bool success = true;

        if (x <= leftVal) {
            success = insert_helper(&(*from)->left, sample, &bub);
        } else {
            success = insert_helper(&(*from)->right, sample, &bub);
        }

        if(!success) {
            
        }

        if (!can_mix((*from)->left->d, (*from)->right->d)) {
            if (bub != nullptr) {
                std::cerr << "variance problem and bubble: " << std::endl;
                std::cerr << (*from)->left->d << " " << (*from)->right->d << " " << bub->d << std::endl;
                std::cerr << mix((*from)->left->d,  (*from)->right->d) << " "
                          << mix((*from)->left->d,  bub->d) << " "
                          << mix((*from)->right->d, bub->d) << std::endl;

                print_helper(std::cerr, 0, (*from)->left);
                print_helper(std::cerr, 0, (*from)->right);
                print_helper(std::cerr, 0, bub);
                
                throw std::logic_error("variance problem and bubble!!!");
            } 

            node * n = *from;
            if ((*from)->left->d.population_var() > (*from)->right->d.population_var()) {
                *from = n->right;
                *bubble = n->left;
            } else {
                *from = n->left;
                *bubble = n->right;
            }
            delete n;
            return true;
        } 
        
        (*from)->d = mix((*from)->left->d, (*from)->right->d);

        // if we hae nothing to bubble, we are done!
        if (bub == nullptr) {
            return true;
        }

        // ok we have a bubble to deal with
        if (can_mix((*from)->d, bub->d)) {
            *from = new node(*from, bub);
            return true;
        }

        // bubble it up again!
        *bubble = bub;
        return true;
    }


    // void insert_helper(node ** from, double sample, node ** bubble) {
    //     if (*from == nullptr) {
    //         *from = new node(sample);
    //         return;
    //     }

    //     if ((*from)->left == nullptr) {
    //         node * r = new node(sample); // r->left/right == nullptr
    //         node * l = *from;            // l->left/right == nullptr == *from
    //         node * n = nullptr;          // n->left == l, n->right == r (dependent on coin flip)
    //         // flip a coin to decide which node to make right versus left
    //         if (uniform(gen) <= 0.5) {
    //             n = new node(r, l);
    //         } else {
    //             n = new node(l, r);
    //         }
    //         // splice in the new node
    //         *from = n;                   // parent node (either right or left) == n
    //         return;
    //     } 


    //     // calculate the probability based on the density
    //     double leftVal = (*from)->left->d.density(sample);
    //     double rightVal = (*from)->right->d.density(sample);

    //     // flip a weighted coin
    //     double x = (leftVal + rightVal) * uniform(gen);
    //     bool insertLeft = (x <= leftVal);

    //     std::cout << "can mix: " << (*from)->left->d << " " << (*from)->right->d << " -> " << can_mix((*from)->left->d, (*from)->right->d) << std::endl;

    //     // std::cout << "PRE left->var: " << vl << std::endl;   

    //     node * bub = nullptr;

    //     if (insertLeft) {
    //         insert_helper(&(*from)->left, sample, &bub);
    //     } else {
    //         insert_helper(&(*from)->right, sample, &bub);
    //     }

    //     if (bub == nullptr) {
    //         // either we can mix or we can't
    //         if (can_mix((*from)->left->d, (*from)->right->d)) {
    //             (*from)->d = mix((*from)->left->d, (*from)->right->d);
    //         } else {
    //             if ((*from)->left->d.var() > (*from)->right->d.var()) {
    //                 // pick the greater variance and bubble it up
    //                 *bubble = (*from)->left;
    //                 node * n = *from;
    //                 *from = (*from)->right;
                    
    //                 delete n;
    //             } else {
    //                 *bubble = (*from)->right;
    //                 node * n = *from;
    //                 *from = (*from)->left;

    //                 delete n;
    //             }
    //         }
    //     } else {
    //         // we need to decide between three nodes, left, right and bub
    //         // sort the three by variance

    //         node * n[] = { (*from)->left, (*from)->right, bub };
    //         std::sort(n, n + 3, [](node * a, node * b) -> bool { return a->d.var() > b->d.var(); });

    //         bool can01 = can_mix(n[0]->d, n[1]->d);
    //         bool can02 = can_mix(n[0]->d, n[2]->d);
    //         bool can12 = can_mix(n[1]->d, n[2]->d);

    //         // is it possible that we can't mix any of them?!
    //         if (!can01 && !can02 && !can12) {
    //             std::cout << "nothing can mix: \n" << n[0]->d << "\n" << n[1]->d << "\n" << n[2]->d << std::endl;
    //             std::cout << "mixes: \n01 - " << mix(n[0]->d, n[1]->d) << "\n"
    //                       << "02 - " << mix(n[0]->d, n[2]->d) << "\n" 
    //                       << "12 - " << mix(n[1]->d, n[2]->d) << std::endl;

    //             print_helper(std::cout, 0, n[0]);
    //             print_helper(std::cout, 0, n[1]);
    //             print_helper(std::cout, 0, n[2]);

    //             throw std::logic_error("nothing can mix, AHHHHH!");
    //         }
    //     }


    //     if ((*from)->d.var() < (*from)->left->d.var()) std::cout << "left variance greater" << std::endl;
    //     if ((*from)->d.var() < (*from)->right->d.var()) std::cout << "right variance greater" << std::endl;
    // }

    void print_helper(std::ostream & os, int depth, node * n) const {
        if (n == nullptr) return;

        for(int d = 0; d < depth; d++) os << "  ";
        os << "[ " << n->d.mean << " ; " << n->d.population_std() << " # " << n->d.count << " -- " << n->error2() << " ]" << std::endl;

        print_helper(os, depth + 1, n->left);
        print_helper(os, depth + 1, n->right);
    }

    void delete_helper(node * from) {
        if (from == nullptr) return;

        delete_helper(from->left);
        from->left = nullptr;
        delete_helper(from->right);
        from->right = nullptr;
        delete from;
    }

public:
    cluster() 
        : root(nullptr), uniform(0,1) 
    { }

    ~cluster() {
        delete_helper(root);
    }

    void insert(double sample) {
        node * bub = nullptr;
        node * s = new node(sample);

        std::set<node*> inserted;


        for(int i = 0; true; i++) {
            std::cout << "try #" << i << std::endl;
            insert_helper2(&root, s, &bub);
            inserted.insert(s);

            // bubbled all the way up-- we need to add a new root
            if (bub != nullptr) {
                if (bub->d.population_var() > root->d.population_var()) {
                    std::swap(root, bub);
                }

                if (inserted.find(bub) != inserted.end()) {
                    throw std::logic_error("trying to insert the same node twice!");
                }
                
                s = bub;
                bub = nullptr;
            } else {
                break;
            }
        }
    }
};

std::ostream & operator<<(std::ostream & os, cluster const & c) {
    c.print_helper(os, 0, c.root);
    return os;
}

int main(int ac, char ** av) {
    dist d0(0, 1, 10);
    dist d1(0, 1, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;

    d0 = dist(0, 2, 10);
    d1 = dist(0, 2, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;
  
    d0 = dist(0, 3, 10);
    d1 = dist(0, 3, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;

    d0 = dist(1, 1, 10);
    d1 = dist(1, 1, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;

    d0 = dist(-5, 5, 10);
    d1 = dist(-5, 5, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;

    d0 = dist(1, 0.5, 10);
    d1 = dist(1, 0.5, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;

    d0 = dist(-1, 1, 10);
    d1 = dist(1, 1, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;

    d0 = dist(-2, 1, 10);
    d1 = dist(2, 1, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;

    d0 = dist(0, 1, 10);
    d1 = dist(0, 2, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;

    d0 = dist(0, 1, 10);
    d1 = dist(0, 0.5, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;


    d0 = dist(-1e100, 1, 10);
    d1 = dist(1e100, 1, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;


    d0 = dist(-1e100, 2, 10);
    d1 = dist(1e100, 2, 10);

    std::cout << d0 << " mix " << d1 << "->" << mix(d0, d1) << " -- " << comp_dist(d0, d1) << std::endl;



    cluster c;

    std::minstd_rand gen;
    std::normal_distribution<double> n1(3, 1);
    std::normal_distribution<double> n2(-3, 1);
    std::uniform_int_distribution<int> coin(0, 1);

    for(int x = 0; x < 20; x++) {
        if (coin(gen) == 0) {
            c.insert3(n1(gen));
        } else {
            c.insert3(n2(gen));
        }
    }

    std::cout << c;

    return 0;
}