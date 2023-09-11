//
// Created by Xandru Mifsud on 28/09/2021.
// Largely based on the following blog post by Martin Ankerl (Published: 08/12/2018, Accessed: 28/09/2021)
// https://martin.ankerl.com/2018/12/08/fast-random-bool/
//

#ifndef ALLLSATISFIABILITYSOLVER_RANDOMBOOLGENERATOR_H
#define ALLLSATISFIABILITYSOLVER_RANDOMBOOLGENERATOR_H

#define UNLIKELY(x) __builtin_expect((x), 0)

#include "random"

typedef unsigned long long ull;

/* Based on the implementation of Martin Ankerl, given at https://martin.ankerl.com/2018/12/08/fast-random-bool/.
 *
 * This class defines an efficient templated random boolean generator, where the template typename E corresponds to a
 * (pseudo)-random number generator engine provided by the <random> header; for further details on these, kindly see the
 * following reference: https://en.cppreference.com/w/cpp/numeric/random.
 *
 * While the <random> header provides a number of random integer generators from which we can stochastically generate 0
 * or 1, these are largely inefficient since they are designed to cater for more general scenarios. Hence the need for
 * k This is especially important since in the worst
 * case of the Algorithmic Lovasz Local Lemma, the number of random boolean re-samples is exponential in the number of
 * variables - hence we require every performance boost possible.
 */
template <typename E>
class RBG{
    public:
        explicit RBG (E& engine) {
            this->engine = engine;
        }

        bool sample() {
            // In case m_rand is equal to initial seed, re-sample (inefficiently) using uniform_int_distribution.
            if (UNLIKELY(1 == m_rand)) {
                m_rand = std::uniform_int_distribution<ull>{}(engine) | s_mask_left1;
            }

            bool const ret = m_rand & 1;
            m_rand >>= 1;
            return ret;
        }

    private:
        const ull s_mask_left1 = ull(1) << (sizeof(ull) * 8 - 1);
        ull m_rand = 1;
        E engine;
};

#endif //ALLLSATISFIABILITYSOLVER_RANDOMBOOLGENERATOR_H
