#include <iostream>
#include <cmath>
#include <random>
#include <limits>
#include <cmath>
#include <utility>
#include <algorithm>

#define SQRT2PI (2.506628274631001)
#define SQRTPI (1.772453851)
#define SQRT2 (1.414213562373095)

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

    double density(double sample) const {
        // if we have a single sample, the density everywhere is 1/sqrt(2 * PI) because it has infinite variance
        if (count == 1) {
            return 1 / SQRT2PI;
        }

        double delta = sample - mean;
        // should we include the count or not?
        // it's going to favor larger clusters, so maybe we should ignore it
        return std::exp(-delta * delta / sample_var() / 2.) / SQRT2PI / sample_std();
    }

    dist() : mean(0), m2(0), count(0) { }
    dist(double mean) : mean(mean), m2(0), count(1) { }
    dist(double mean, double stddev, int count) : mean(mean), m2(stddev * stddev * count), count(count) { }

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

double comp_dist(dist const & left, dist const & right) {
    return gaussian(left.mean, left.population_var() + right.population_var(), right.mean);
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
        { }

        node(node * left, node * right, dist d) 
          : left(left), right(right), d(d)
        { }

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
        insert_helper(&root, s, &bub);

        // bubbled all the way up-- we need to add a new root
        if (bub != nullptr) {
            node * n = new node(bub, root);
            root = n;
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



    cluster c;

    std::minstd_rand gen;
    std::normal_distribution<double> n1(3, 1);
    std::normal_distribution<double> n2(-3, 1);
    std::uniform_int_distribution<int> coin(0, 1);

    for(int x = 0; x < 20; x++) {
        if (coin(gen) == 0) {
            c.insert(n1(gen));
        } else {
            c.insert(n2(gen));
        }
    }

    std::cout << c;

    return 0;
}