#include "part_init.h"
Particle part_init(int N){
int nn;
Particle pHead = (Particle)malloc(sizeof(particle));
Particle pTail = pHead;
pTail->pNext = NULL;
for(nn=0;nn<N;nn++){
Particle pNew = (Particle)malloc(sizeof(particle));
pTail->pNext = pNew;
pNew->pNext = NULL;
pTail = pNew;
}
pHead->ID = -1;
return pHead;
}


Particle part_add(Particle pHead, Particle pAdd){
pAdd->pNext=pHead;
pHead=pAdd;
return pHead;
}

Particle part_del(Particle pHead, Particle pDel){
int k=0;
Particle pH, pM;
pH=pHead;
while((pH!=pDel)||pH){
pM=pH;
pH=pH->pNext;
k=1;
}
if(k&&pH){pM->pNext=pH->pNext;}
else if(!k){pHead=pHead->pNext;}
return pHead;
}