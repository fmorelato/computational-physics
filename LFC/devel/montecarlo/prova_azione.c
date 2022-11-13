#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../include/azione.h"
#include "../../include/global.h"


int main() {

float s;
int i;
    for(i=0; i<N; i++){
        scanf("%lf",&xx[i]);
    }
s = azione();
printf("il valore dell'azione è %f", s );

return 0;

}