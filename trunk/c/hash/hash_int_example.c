#include <stdio.h>            // MUST follow these 5 steps in about the right order !!!
#include "khash.h"            // 1. include khash.h

typedef int key_type;         // 2. typde key_type
typedef char * val_type;      //    and val_type 
#define MISSING_VALUE "NULL"  // 3. define MISSING_VALUE. to match val_type
HASH_INIT_INT_KEYS(val_type)  // 4. initialize. (either HASH_INIT_STR_KEYS or HASH_INIT_INT_KEYS)

#include "khash_sugar.h"      // 5. include khash_sugar.h


int main() {
    hash *h = hash_new();
    int key = 1234334342;
    char *value = "3.11";
    hash_put(h, key, value);

    printf("%s\n", hash_get(h, key));
    
    printf("%s\n", hash_get(h, 123));
    printf("%s\n", hash_get(h, key));
    printf("%i\n", hash_contains(h, key));
    printf("%i\n", hash_contains(h, 22));

    char *v2 = "hello";
    hash_put(h, 33, v2);

    printf("%i\n", hash_size(h));
    printf("%s\n", hash_get(h, key));
    hash_delete(h, key);
    printf("%s\n", hash_get(h, key));

    hash_clear(h);
    hash_free(h);
    return 0;
}

