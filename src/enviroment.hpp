
#define TO_STRING2(X) #X
#define TO_STRING(X) TO_STRING2(X)
#define INSTALLATION_PATH TO_STRING(BIN_PATH)

/* Test for GCC > 4.9.4 */
#define GCC_VERSION (__GNUG__ * 10000 \
                     + __GNUG_MINOR__ * 100 \
                     + __GNUG_PATCHLEVEL__)
#if GCC_VERSION > 40904
    /* std::valarrays do not work well with gcc-4.9.5 and older. When
     * using the following constructor
     * std::valarray<float> arr(size, 0);
     * there is a SEGFAULT when doing arr[1], with 1 < size.
     * The workaround is to perform an explicit resize:
     * std::valarray<float> arr;
     * arr.resize(size);
     * Also, the std::begin and std::end methods are not specialized
     * for valarray according to gcc-4.9.5. More problems arise when
     * subsetting an array using a slice and trying to sum it up.
     */
    #define __GOOD_VALARRAYS
#endif