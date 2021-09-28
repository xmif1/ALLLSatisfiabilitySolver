//
// Created by Xandru Mifsud on 28/09/2021.
// Largely based on the following blog post by Martin Ankerl (Published: 08/12/2018, Accessed: 28/09/2021)
// https://martin.ankerl.com/2018/12/08/fast-random-bool/
//

#ifndef ALLLSATISFIABILITYSOLVER_RANDOMBOOLGENERATOR_H
#define ALLLSATISFIABILITYSOLVER_RANDOMBOOLGENERATOR_H

#define UNLIKELY(x) __builtin_expect((x), 0)

template <typename E, typename U = uint64_t> class RBG{
    public:
        explicit RBG(E& engine){
            this->engine = engine;
        }

        bool sample(){
            if(UNLIKELY(1 == m_rand)){
                m_rand = std::uniform_int_distribution<U>{}(engine) | s_mask_left1;
            }

            bool const ret = m_rand & 1;
            m_rand >>= 1;
            return ret;
        }

    private:
        static constexpr const U s_mask_left1 = U(1) << (sizeof(U) * 8 - 1);
        U m_rand = 1;
        E engine;
};

#endif //ALLLSATISFIABILITYSOLVER_RANDOMBOOLGENERATOR_H
