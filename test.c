#include <stdio.h>
#include <stdlib.h>

int main(){

    FILE *fPtr;
    char s[50];

    fPtr = fopen("oldname.txt", "r");
    
    fgets(s, 50, fPtr);
    printf(s,s[0]);
    fgets(s, 50, fPtr);
    printf(s);
    
    printf("\n");
 
    fclose(fPtr);
     
    return 0;


}
