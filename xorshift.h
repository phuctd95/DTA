#ifndef XORSHIFT_H
#define XORSHIFT_H

#include <stdint.h>
#include <vector>

class xorshift
{
    private:
        uint64_t s[16];
        int32_t p;
        uint64_t get(void);

    public:
        xorshift();
        virtual ~xorshift();
        void init_seed(int32_t p_init, std::vector <uint64_t> &s_init);
        uint64_t get_max( uint64_t m);
        double get_real();
};

#endif // XORSHIFT_H
