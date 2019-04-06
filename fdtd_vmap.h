#ifndef _FDTD_VMAP_H
#define _FDTD_VMAP_H
#include "math.h"
#include "fdtd_grid.h"

/* Map variables a easy form, make a series by yourself */

/* 2d grid */
#define HxG(G, M, N) G->hx[(M) * (SizeYG(G)-1) + (N)]
#define HxoldG(G, M, N) G->hxold[(M) * (SizeYG(G)-1) + (N)]
#define ChxhG(G, M, N) G->chxh[(M) * (SizeYG(G)-1) + (N)]
#define ChxeG(G, M, N) G->chxe[(M) * (SizeYG(G)-1) + (N)]
#define mxG(G, M, N) G->mx[(M) * (SizeYG(G)-1) + (N)]
#define peyzG(G, M, N) G->peyz[(M) * (SizeYG(G)-1) + (N)]
#define peyz_bG(G, M) G->peyz_b[(M)]
#define peyz_cG(G, M) G->peyz_c[(M)]
#define peyz_kG(G, M) G->peyz_k[(M)]






#define HyG(G, M, N) G->hy[(M) * SizeYG(G) + (N)]
#define HyoldG(G, M, N) G->hyold[(M) * SizeYG(G) + (N)]
#define ChyhG(G, M, N) G->chyh[(M) * SizeYG(G) + (N)]
#define ChyeG(G, M, N) G->chye[(M) * SizeYG(G) + (N)]
#define myG(G, M, N) G->my[(M) * SizeYG(G) + (N)]
#define pexzG(G, M, N) G->pexz[(M) * SizeYG(G) + (N)]
#define pexz_bG(G, M) G->pexz_b[(M)]
#define pexz_cG(G, M) G->pexz_c[(M)]
#define pexz_kG(G, M) G->pexz_k[(M)]






#define HzG(G, M, N) G->hz[(M) * (SizeYG(G)-1) + (N)]
#define HzoldG(G, M, N) G->hzold[(M) * (SizeYG(G)-1) + (N)]
#define ChzhG(G, M, N) G->chzh[(M) * (SizeYG(G)-1) + (N)]
#define ChzeG(G, M, N) G->chze[(M) * (SizeYG(G)-1) + (N)]
#define mzG(G, M, N) G->mz[(M) * (SizeYG(G)-1) + (N)]
#define pexyG(G, M, N) G->pexy[(M) * (SizeYG(G)-1) + (N)]
#define peyxG(G, M, N) G->peyx[(M) * (SizeYG(G)-1) + (N)]
#define DensityG(G, M, N) G->density[(M) * (SizeYG(G)-1) + (N)]




#define ExG(G, M, N) G->ex[(M) * SizeYG(G) + (N)]
#define CexeG(G, M, N) G->cexe[(M) * SizeYG(G) + (N)]
#define CexhG(G, M, N) G->cexh[(M) * SizeYG(G) + (N)]
#define jxG(G, M, N) G->jx[(M) * SizeYG(G) + (N)]
#define phyzG(G, M, N) G->phyz[(M) * SizeYG(G) + (N)]
#define phyz_bG(G, M) G->phyz_b[(M)]
#define phyz_cG(G, M) G->phyz_c[(M)]
#define phyz_kG(G, M) G->phyz_k[(M)]









#define EyG(G, M, N) G->ey[(M) * (SizeYG(G)-1) + (N)]
#define CeyeG(G, M, N) G->ceye[(M) * (SizeYG(G)-1) + (N)]
#define CeyhG(G, M, N) G->ceyh[(M) * (SizeYG(G)-1) + (N)]
#define jyG(G, M, N) G->jy[(M) * (SizeYG(G)-1) + (N)]
#define phxzG(G, M, N) G->phxz[(M) * (SizeYG(G)-1) + (N)]
#define phxz_bG(G, M) G->phxz_b[(M)]
#define phxz_cG(G, M) G->phxz_c[(M)]
#define phxz_kG(G, M) G->phxz_k[(M)]







#define EzG(G, M, N) G->ez[(M) * SizeYG(G) + (N)]
#define CezeG(G, M, N) G->ceze[(M) * SizeYG(G) + (N)]
#define CezhG(G, M, N) G->cezh[(M) * SizeYG(G) + (N)]
#define jzG(G, M, N) G->jz[(M) * SizeYG(G) + (N)]
#define phxyG(G, M, N) G->phxy[(M) * SizeYG(G) + (N)]
#define phyxG(G, M, N) G->phyx[(M) * SizeYG(G) + (N)]









#define SizeXG(G) G->sizeX
#define SizeYG(G) G->sizeY
#define TimeG(G) G->time
#define MaxTimeG(G) G->maxTime
#define CdtdsG(G) G->cdtds
#define TypeG(G) G->type
#define lxG(G) G->lx
#define lyG(G) G->ly
#define DtG(G) G->Dt
#define lamdG(G) G->lamd
#define ulG(G) G->ul
#define utG(G) G->ut

#define pmlxlG(G) G->pmlxl
#define pmlxrG(G) G->pmlxr
#define pmlyuG(G) G->pmlyu
#define pmlydG(G) G->pmlyd




/* macros that assume the "Grid" is "g" */

/* 2d grid */
#define Hx(M, N) HxG(g, M, N)
#define Hxold(M, N) HxoldG(g, M, N)
#define Chxh(M, N) ChxhG(g, M, N)
#define Chxe(M, N) ChxeG(g, M, N)
#define Mx(M, N) mxG(g, M, N)
#define Peyz(M, N) peyzG(g, M, N)
#define Peyz_b(M) peyz_bG(g, M)
#define Peyz_c(M) peyz_cG(g, M)
#define Peyz_k(M) peyz_kG(g, M)





#define Hy(M, N) HyG(g, M, N)
#define Hyold(M, N) HyoldG(g, M, N)
#define Chyh(M, N) ChyhG(g, M, N)
#define Chye(M, N) ChyeG(g, M, N)
#define My(M, N) myG(g, M, N)
#define Pexz(M, N) pexzG(g, M, N)
#define Pexz_b(M) pexz_bG(g, M)
#define Pexz_c(M) pexz_cG(g, M)
#define Pexz_k(M) pexz_kG(g, M)





#define Hz(M, N) HzG(g, M, N)
#define Hzold(M, N) HzoldG(g, M, N)
#define Chzh(M, N) ChzhG(g, M, N)
#define Chze(M, N) ChzeG(g, M, N)
#define Mz(M, N) mzG(g, M, N)
#define Pexy(M, N) pexyG(g, M, N)
#define Peyx(M, N) peyxG(g, M, N)
#define Density(M, N) DensityG(g, M, N)




#define Ex(M, N) ExG(g, M, N)
#define Cexe(M, N) CexeG(g, M, N)
#define Cexh(M, N) CexhG(g, M, N)
#define Jx(M, N) jxG(g, M, N)
#define Phyz(M, N) phyzG(g, M, N)
#define Phyz_b(M) phyz_bG(g, M)
#define Phyz_c(M) phyz_cG(g, M)
#define Phyz_k(M) phyz_kG(g, M)





#define Ey(M, N) EyG(g, M, N)
#define Ceye(M, N) CeyeG(g, M, N)
#define Ceyh(M, N) CeyhG(g, M, N)
#define Jy(M, N) jyG(g, M, N)
#define Phxz(M, N) phxzG(g, M, N)
#define Phxz_b(M) phxz_bG(g, M)
#define Phxz_c(M) phxz_cG(g, M)
#define Phxz_k(M) phxz_kG(g, M)




#define Ez(M, N) EzG(g, M, N)
#define Ceze(M, N) CezeG(g, M, N)
#define Cezh(M, N) CezhG(g, M, N)
#define Jz(M, N) jzG(g, M, N)
#define Phxy(M, N) phxyG(g, M, N)
#define Phyx(M, N) phyxG(g, M, N)





#define SizeX SizeXG(g)
#define SizeY SizeYG(g)
#define Time TimeG(g)

#define MaxTime MaxTimeG(g)

#define Lx lxG(g)
#define Ly lyG(g)
#define dt DtG(g)
#define Lamd lamdG(g)
#define Ul ulG(g)
#define Ut utG(g)


#define cpmlxl pmlxlG(g)
#define cpmlxr pmlxrG(g)
#define cpmlyu pmlyuG(g)
#define cpmlyd pmlydG(g)



#endif 