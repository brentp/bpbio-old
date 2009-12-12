/*
 * this is a small, simple set of macros to provide a simple 
 * API to khash.h : http://attractivechaos.wordpress.com/2008/09/02/implementing-generic-hash-library-in-c/
 * LICENSE is MIT, same as khash.h
 * see hash_example for useage.
 */

// MUST follow these 5 steps in about the right order !!!
// ----------------------------------------------------//
// 1. include khash.h
// 2. typde key_type
//    and val_type 
// 3. define MISSING_VALUE. to match val_type
// 4. initialize. (either HASH_INIT_STR_KEYS or HASH_INIT_INT_KEYS)
// 5. include khash_sugar.h



#define hash_init(key_type, val_type) KHASH_MAP_INIT_#key_type(, val_type)
static int __ret;
#define hash_contains(h, key) kh_exist((h), kh_get(, h, (key)))
#define hash_size(h) ((h)->size)
#define hash_clear(h) (kh_clear(, h))
#define hash_free(h) kh_destroy(, (h))
#define hash_delete(h, k) kh_del(, h, kh_get(,h, (key)))
#define hash_new() kh_init()
#define hash khash_t()

inline void hash_put(hash *h, key_type key, val_type value){
    khiter_t k = kh_put(,h, key, &__ret);
    kh_value(h, k) = value;
}

inline val_type hash_get(hash *h, key_type key){
    khiter_t k = kh_get(,h, key);
    return k == kh_end(h) ? MISSING_VALUE : kh_value(h, k);
}
