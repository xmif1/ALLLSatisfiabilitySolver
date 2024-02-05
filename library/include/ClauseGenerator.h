//
// Created by Xandru Mifsud on 10/09/2023.
//

#ifndef ALLLSATISFIABILITYSOLVERMAIN_CLAUSEGENERATOR_H
#define ALLLSATISFIABILITYSOLVERMAIN_CLAUSEGENERATOR_H

#include "Clause.h"

using namespace std;

template<class T, class Enable = void>
class ClauseGenerator {}; // Default specialisation

template <class T> // Specialisation to integral types
class ClauseGenerator<T, typename enable_if<is_integral<T>::value>::type> {

public:
    using ClauseArray = typename Clause<T>::ClauseArray;
    typedef unsigned short int t_id_T;

    T n_clauses;

    ClauseGenerator(Clause<T>* (*getEnumeratedClause)(T, t_id_T), t_id_T t_id, T n_clauses, T base_offset, T batch_size) {
        this->getEnumeratedClause = getEnumeratedClause;
        this->t_id = t_id;
        this->n_clauses = n_clauses;
        this->base_offset = base_offset;
        this->batch_size = batch_size;
    }

    ClauseArray* yieldRandomUNSATClauseBatch(const bool* var_arr) {
        if (finished_yielding) {
            reset();
        }

        auto clauses = new ClauseArray();

        T n = batch_size;
        if (n_yielded_clauses + batch_size >= n_clauses) {
            n = n_clauses - n_yielded_clauses;
        }

        for (T i = 0; i < n; i++) {
            c = (c + P) % n_clauses;
            auto clause = getEnumeratedClause(base_offset + c, t_id);

            if (clause != nullptr) {
                if (clause->is_not_satisfied(var_arr)) {
                    clauses->push_back(clause);
                } else { // Memory management
                    clause->literals->clear();
                    delete clause->literals;
                    delete clause;
                }

                n_yielded_clauses++;
            } else {
                cerr << "WARNING: Clause generator went out of range and yielded nullptr." << endl;

                finished_yielding = true;
                break;
            }
        }

        if (n_yielded_clauses == n_clauses) {
            finished_yielding = true;
        }

        return clauses;
    }

    Clause<T>* yieldNextClause() {
        if (finished_yielding) {
            reset();
        }

        auto clause = getEnumeratedClause(base_offset + n_yielded_clauses, t_id);

        if (clause != nullptr) {
            n_yielded_clauses++;
            if (n_yielded_clauses == n_clauses) {
                finished_yielding = true;
            }

            return clause;
        } else {
            cerr << "WARNING: Clause generator went out of range and yielded nullptr." << endl;

            finished_yielding = true;
            return nullptr;
        }
    }

    bool has_finished_yielding() {
        return finished_yielding;
    }

    void reset() {
        n_yielded_clauses = 0;
        finished_yielding = false;
    }

    private:
        Clause<T>* (*getEnumeratedClause)(T, unsigned short int);
        T batch_size, base_offset;
        t_id_T t_id{};

        const uint64_t P = 9223372036854775783;
        T c = 0;

        bool finished_yielding = false;
        T n_yielded_clauses = 0;
};

#endif //ALLLSATISFIABILITYSOLVERMAIN_CLAUSEGENERATOR_H
