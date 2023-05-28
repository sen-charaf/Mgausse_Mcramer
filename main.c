#include <stdio.h>
#include <stdlib.h>
#include <time.h>
typedef struct{
int nmr; //Numérateur
int dnm; //Dénominateur
}rational;
rational add_r(rational a,rational b); //addition
rational mult_r(rational a,rational b);//multiplication
rational div_r(rational a,rational b);//division
rational redc_r(rational a);    //réductoin

void menu_0(); //Le premier interface qui s'affiche
int menu_1(); //La deusieme interface qui s'affiche
int menu_2(); //La troisieme interface qui s'affiche
rational** allo_tab_2D(const int nbr); //L'allocation dynamique d'une tab 2D
rational* allo_pivo(const int nbr); //L'allocation dynamique d'une tab 1D
void Init_tab_2D(rational** T,const int nbr,const int rep); //L'initialisation de la tab 2D
void Aff_tab_2D(const rational** T,const int nbr); //L'affichage de la tab 2D
void Init_pivo(rational* R,const rational** T,const int nbr);//L'initialisation de la tab 1D
int Modi_tab_2D(rational** T,const int nbr); /*Affecter des modifications sur
                                            les elements de la tab 2D*/
int Switch_tab_2D(rational** T,const int nbr); /*switchi les linges
                                                                    de la tab 2d*/
int Deter_tab_2D(const rational* R,const int nbr,rational pro); /*Calculer la déterminant
                                                            avec la methode de Gauss*/
rational Cramer_det(const rational** T,const int nbr,int k,int j);/*Calculer la déterminant
                                                            avec la methode de Cramer*/
void Init_cramer(rational** R,const rational** T,int l,int m,const int nbr);
rational** deallo_tab_2D(rational** T,const int nbr); //La dealloca dynamique de la tab 2D
rational* deallo_pivo(rational* R); //La deallocation dynamique de la tab 1D
/*
Pour bien comprendre les fonctions au-dessus vous pouvez aller
 à leur définition au dessou de la fonction main
 (ils sont en le même ordre de la déclaration)
*/

int main()
{
    /*Block de déclaration:*/
    /*====================*/
    int nb,i,j,k,opt,ind0,mthd;   //opt:option  //nb:nombre //ind0:l'indice d'annulation de pivot
    char rep;               // rep:reponse
    rational tmp1,tmp2,det,detc,ops;  //tmp:temporelle //det:déterminant
    //ops:opposée
    rational** M;             // La matrice
    rational* Piv;            // Le pivot (un tableau 1D)
    /*====================*/


    /*cette bouqule signifie
    que ce programme va se
    relancer tant que l'utilisateur
    veut continuer l'exécution*/
    /*====================*/
    do
    {
        ops.nmr=-1;
        det.nmr=det.dnm=ops.dnm=1;
        ind0=1;
        menu_0();

        /*
        opt va contenir le choix saisit par l'utilisateur
        pour determiner la manière d'initialisation
        de la taille de la matrice (1) ,et la manière de
        s'initialisation (2)
        */
        mthd=menu_1();
        opt=menu_2();
        if(opt==1||opt==2)  // (1)
        {
            do
            {
                printf("\t\t\tDonner la taille de la matrice: ");
                scanf("%d",&nb);
            }
            while(nb<=0);
        }
        else
            nb=4;
        if(mthd==2){
            M=allo_tab_2D(nb);
            Init_tab_2D(M,nb,opt);
            system("cls");
            Aff_tab_2D(M,nb);
            Modi_tab_2D(M,nb);
            Aff_tab_2D(M,nb);
            printf("\nDeterminant= ");
            detc=Cramer_det(M,nb,0,0);
            if(detc.dnm==1||detc.nmr==0)
                printf(" %d",detc.nmr);
            else {
                detc=redc_r(detc);
                printf(" %d/%d",detc.nmr,detc.dnm);
            }
        }else{
        M=allo_tab_2D(nb);
        Piv=allo_pivo(nb);
        Init_tab_2D(M,nb,opt); // (2)
        system("cls");
        Aff_tab_2D(M,nb);
        Modi_tab_2D(M,nb);
        Aff_tab_2D(M,nb);
        /*Ce linge permer de mulutplier la déterminant
        par -1 à chaque réalisation de switch*/
        det.nmr=Switch_tab_2D(M,nb);
        Init_pivo(Piv,M,nb);
        for(j=0; j<nb; j++)
        {
            /*initialiser i avec j+1 pour avoir
            une matrice triangulait et pour
            ne pas annuler les elements du diagonal*/
            for(i=j+1; i<nb; i++)
            {
                if(!M[i][j].nmr)
                    continue; //ignorer les cas où il existe un zéro
                tmp1=M[i][j];  //mémoriser l'élément au cours d'élimination
                for(k=j; k<nb; k++)
                {
                    M[i][k]=add_r(M[i][k],mult_r(div_r(mult_r(ops,tmp1),Piv[j]),M[j][k]));//formule génerale de l'élimination
                    M[i][k]=redc_r(M[i][k]);
                    if(i==k)
                        Piv[i]=M[i][k]; //Updating le pivot
                    if(Piv[i].nmr==0 && i==0)
                    {
                        /*Si le pivot était nul par élimination
                        on va sortir des boucles  et le déterminant est égal à zéro*/
                        ind0=i; //ind0 pour mémoriser i et l'utiliser dans (3) et (4)
                        break;
                    }
                }
                if(Piv[ind0].nmr==0 && ind0==j)  //(3)
                {
                    Aff_tab_2D(M,nb);

                    break;
                }
                Aff_tab_2D(M,nb);
                tmp2=div_r(mult_r(ops,tmp1),Piv[j]);
                if(tmp2.dnm==1||tmp2.nmr==0)
                    printf("\tL%d <-- L%d + %d * L%d",i+1,i+1,tmp2.nmr,j+1);
                else
                    printf("\tL%d <-- L%d + %d/%d * L%d",i+1,i+1,tmp2.nmr,tmp2.dnm,j+1);
                getch();
            }
            if(Piv[ind0].nmr==0 && ind0==j) // (4)
                break;
        }
        Deter_tab_2D(Piv,nb,det);
        }
        printf("\n");
        getch();
        system("cls");

        /*Après la deallocation de M et Piv on l'affecter
        le pointeur null pour éviter les erreurs*/
        M=deallo_tab_2D(M,nb);
        if(mthd==1)
            Piv=deallo_pivo(Piv);
        /*Domander l'autorisation de quitter le programme
         ou bien de le relancer une autrefois*/
        printf("\t\t\tVoullez vous Retester Le programme?\n\t\t\t\t\t O: Oui / N: Non\n\n");
        do
        {
            printf("\t\t\t");
            scanf(" %c",&rep);
        }
        while(rep!='O'&&rep!='o'&&rep!='N'&&rep!='n');
        system("cls");
    }
    while(rep=='O'||rep=='o');
    /*====================*/

    return 0;
}
rational add_r(rational a,rational b)
{
    rational tmp;
    if(a.dnm==b.dnm){
        /*On n'a pas besoin d'uniter le dnm*/
        tmp.nmr=a.nmr+b.nmr;
        return tmp;
    }else{
        /*On a besoin d'uniter le dnm*/
        /*On va muliplier le nmr de "a" avec
        le dmn de "b"et ainsi de suite*/
        tmp.nmr=a.nmr;
        tmp.dnm=a.dnm;
        a.nmr*=b.dnm;
        a.dnm*=b.dnm;
        b.nmr*=tmp.dnm;
        b.dnm*=tmp.dnm;
        tmp.nmr=a.nmr+b.nmr;
        tmp.dnm=a.dnm;
        return tmp;
    }
}
rational mult_r(rational a,rational b)
{
    /*Muliplier le nmr de "a" avec le nmr de
    "b" et ainsi de suite pour le dnm*/
    rational tmp;
    tmp.nmr=a.nmr*b.nmr;
    tmp.dnm=a.dnm*b.dnm;
    return tmp;
}
rational div_r(rational a,rational b)
{
    /*Muliplier le nmr de "a" avec le dnm de
    "b" et le dnm de "a" avec le nmr de "b"*/
    rational tmp;
    tmp.nmr=a.nmr*b.dnm;
    tmp.dnm=a.dnm*b.nmr;
    return tmp;
}
rational redc_r(rational a)
{
    /*Cette fonction est permet de réduire un nombre rationelle
    selon les cas au-dessous*/
    rational tmp;
    tmp.nmr=tmp.dnm=1;
    if(a.nmr==a.dnm)
        /*Si nmr==dnm la fonctions reteurne 1*/
        return tmp;
    if(a.nmr<0 && a.dnm<0)
    {
        /*On élémine les singes moins
        si ils sont au le nmr et le dnm*/
        a.nmr=abs(a.nmr);
        a.dnm=abs(a.dnm);
    }

    int i=2,PGCD=1,max;
    /*Determiner le PGCD"Plus grand diviseur commun"
    et diviser le nmr et le dnm par lui*/
    if(abs(a.nmr)>abs(a.dnm))
        //On va prendre le maximume rang
        max=a.nmr;
    else
        max=a.dnm;
    do{
        if((a.nmr%i)==0 && (a.dnm%i)==0)
            PGCD=i;
        i++;
    }while(i<=abs((max)/2));
    tmp.nmr=a.nmr/PGCD;
    tmp.dnm=a.dnm/PGCD;
    return tmp;
}
void menu_0()
{
    /* Cette fonction est simplement permet d'afficher un premier menu de
    programme et de donner la choisit de continuer ou exiter le programme*/
    char rep;
    system("cls");
    printf("\n\t\t======================================================");
    printf("\n\t\tBien venue a ce programme qu'il est permer de calculer\n\t\tLa determinant d'une matrice Carree d'ordre quelconque\n\t\tAvec la methode de Gauss.");
    printf("\n\t\t======================================================");
    printf("\n\n\t\t\tRealises par: Charaf-eddine Kaouri");
    printf("\n\n\t\t\t      De 11/11/22 A 12/11/22");
    printf("\n\t\t------------------------------------------------------");
    printf("\n\n\t\tCe programme est capable de vous donner la possibiliter:");
    printf("\n\n\t\t-Initaliser une matrice d'une maniere:");
    printf("\n\n\t\t\t.Random (avec determiner le rang de rand)\n\n\t\t\t.Manuelle\n\n\t\t\t.D'importer un exemple typique (A4x4)");
    printf("\n\n\t\t-Changer les elements de matrice avant le calule");
    printf("\n\n\t\t-Calculer et Afficher la determinant de la matrice");
    printf("\n\n\t\t\tVoullez vous continuer?");
    printf("\n\n\t\t\tO: Oui / N: Non");
    do
    {
        printf("\n\t\t\t");
        scanf(" %c",&rep);
    }
    while(rep!='O'&&rep!='o'&&rep!='N'&&rep!='n');
    if(rep=='N'||rep=='n')
        exit(5);
    system("cls");
}
int menu_1()
{
    /*Cette fonction est permet de demander à l'utilisateur de choisir
     une methode (Gauss,Cramer)
    Et retourner un variable "rep" qu'il contient la réponse*/
    int rep;
    printf("\n\n\t\t\t\tMETHODE:");
    printf("\n\t--------------------------------------------------------------------------");
    printf("\n\n\tPour choisir l'une des methodes suivantes entree le chiffre convenable: ");
    printf("\n\n\t\tGauss: 1\n\n\t\tCramer: 2");
    do
    {
        printf("\n\n\t\t");
        scanf("%d",&rep);
    }
    while(rep!=1&&rep!=2);
    system("cls");
    return rep;
}
int menu_2()
{
    /*Cette fonction est permet d'afficher un menu contiene  plusieurs
    choix d'initialisation, demander à l'utilisateur de choisir une.
    Et retourner un variable "rep" qu'il contient la réponse*/
    int rep;
    system("cls");
    printf("\n\n\t\t\t\tINITIALISATION:");
    printf("\n\t--------------------------------------------------------------------------");
    printf("\n\n\tPour choisir l'une des initialisations suivantes entree le chiffre convenable: ");
    printf("\n\n\t\tRandom: 1\n\n\t\tManuelle: 2\n\n\t\tun exemple typique (A4x4): 3");
    do
    {
        printf("\n\n\t\t");
        scanf("%d",&rep);
    }
    while(rep!=1&&rep!=2&&rep!=3);
    system("cls");
    return rep;
}
rational** allo_tab_2D(const int nbr)
{
    /*Cette fonction est permet d'allouer et de calculer
    un espace de mémoire pour un tableau de 2D et de
    retourner un pointeur si elle reausit sinon
    elle va quitter le programme */
    int i;
    rational** T;
    T=(rational**)malloc(nbr*sizeof(rational*));
    if(T==NULL)
    {
        printf("\nERROR");
        exit(1);
    }
    for(i=0; i<nbr; i++)
    {
        T[i]=(rational*)malloc(nbr*sizeof(rational));
        if(T[i]==NULL)
        {
            printf("\nERROR");
            exit(2);
        }
    }
    /*Le calcule de l'espace allouée*/
    //printf("Matrice:\n%zu Bytes  Successfully allocated!!\n",nbr*sizeof(rational*)+nbr*sizeof(rational));
    return T;
}
rational* allo_pivo(const int nbr)
{
    /*Cette fonction est permet d'allouer et de calculer
    un espace de mémoire pour un tableau de 1D et de
    retourner un pointeur si elle reausit sinon
    elle va quitter le programme*/
    rational* R;
    R=(rational*)malloc(sizeof(rational)*nbr);
    if(R==NULL)
    {
        printf("\nERROR");
        exit(3);
    }
    /*Le calcule de l'espace allouée*/
    //printf("Pivot:\n%zu Bytes  Successfully allocated!!\n",nbr*sizeof(rational));
    system("cls");
    return R;
}
void Init_tab_2D(rational** T,const int nbr,const int rep)
{
    /* Cette fonction est permet d'initialiser
    la matrice "T" par des maniers différents selon
    l'argument "rep"*/
    int i,j,rmax,rmin;
    // rmax/rmin:la valeur maximale/minimale du nombre randome
    int Y[4][4]= {{1,0,3,2},{0,1,-1,0},{2,2,1,1},{-1,0,1,3}};
    rational R[4][4];
    for(i=0; i<4; i++)
            for(j=0; j<4; j++){
                R[i][j].nmr=Y[i][j];
                //R[i][j].dnm=1;
            }
    //declaration de l'exemple typique
    for(i=0; i<nbr; i++)
            for(j=0; j<nbr; j++)
                T[i][j].dnm=1;
    switch (rep)
    {
    case 1:
        srand(time(NULL)); //initialisation de seed de rand chaque segand
        printf("\n\n\t\t");
        printf("Donner le rang de les nombres randoms\n\n\t\tLe minimume: ");
        scanf("%d",&rmin);
        printf("\n\n\t\t");
        printf("Le maximume: ");
        scanf("%d",&rmax);
        for(i=0; i<nbr; i++)
            for(j=0; j<nbr; j++)
                T[i][j].nmr=(rand()%(rmax-rmin+1))+rmin;
        break;
    case 2:
        for(i=0; i<nbr; i++)
            for(j=0; j<nbr; j++)
            {
                printf("Donner l'elemant M[%d][%d]= ",i,j);
                scanf("%d",&T[i][j].nmr);
            }
        break;
    case 3:
        for(i=0; i<nbr; i++)
            for(j=0; j<nbr; j++)
                T[i][j].nmr=R[i][j].nmr;
        break;
    }

}
void Aff_tab_2D(const rational** T,const int nbr)
{
    /* Cette fonction est permet d'afficher la matrice */
    int i,j;
    printf("\n");
    printf("\n");
    for(i=0; i<nbr; i++)
    {
        printf("\n| ");
        for(j=0; j<nbr; j++){
            if(T[i][j].dnm==1||T[i][j].nmr==0)
                printf(" %3d ",T[i][j].nmr);
            else
                printf(" %3d/%d ",T[i][j].nmr,T[i][j].dnm);
        }
        printf("\t| ");
    }
}
void Init_pivo(rational* R,const rational** T,const int nbr)
{
    /* Cette fonction est permet d'initialiser le pivot
    (un tableau de 1D) qu'il contient les elements du diagonal */
    int i,j;
    for(i=0; i<nbr; i++)
        for(j=0; j<nbr; j++)
            if(i==j){
                R[i].nmr=T[i][j].nmr;
                R[i].dnm=T[i][j].dnm;
            }
}
int Modi_tab_2D(rational** T,const int nbr)
{
    /*Cette fonction est permet de demander à l'utilisateur de choisir
     S'il veut modifier la matrice, s'il répond par oui la fonction
     va demander à l'utilisateur de choisir la position de l'element
    qu'il veut modifier et leur nouvelle valeur*/
    char rep;
    int i,j;
    rational nvl; //nouvelle valeur
    printf("\n\n\t\t\tVoullez vous modifier un elements?");
    printf("\n\n\t\t\tO: Oui / N: Non");
    printf("\n\n\t\t");
    do
    {
        printf("\n\n\t\t");
        scanf(" %c",&rep);
    }
    while(rep!='O'&&rep!='o'&&rep!='N'&&rep!='n');
    if(rep=='N'||rep=='n')
    {
        system("cls");
        return 0;
    }
    do
    {
        do
        {
            printf("Entrer i de M[i][j]: ");
            scanf("%d",&i);
        }
        while(i>=nbr);
        if(i<0)
            break;
        do
        {
            printf("Entrer j de M[%d][j]: ",i);
            scanf("%d",&j);
        }
        while(j>=nbr);
        if(j<0)
            break;
        printf("Donner la nouvelle valeur de M[%d][%d]= ",i,j);
        scanf("%d",&nvl.nmr);
        T[i][j].nmr=nvl.nmr;
    }
    while(1);
    /*la fonction n'arrête pas tant que l'utilisateur entre
    des valeurs positives*/
    system("cls");
    return 0;
}
int Switch_tab_2D(rational** T,const int nbr)
{
    //Permer de "switcher" les linges pour eviter le 0 dans le diagonal (au cas de l'existence)
    rational Conv[nbr];            // Conv:converteur
    int i,j,k,rep,etat,tmp;
    rep=1;
    for(i=0; i<nbr-1; i++)
        for(j=0; j<nbr-1; j++){
    if (i==j && !T[i][j].nmr)
    {

            tmp=0;
            while(!T[tmp][j].nmr){
                if(tmp>=nbr-1)
                return 1;
                tmp++;}
            for(k=0; k<nbr; k++)
            {
                Conv[k].nmr=T[i][k].nmr;
                T[i][k].nmr=T[tmp][k].nmr;
                T[tmp][k].nmr=Conv[k].nmr;
                etat=tmp;
            }


    rep*=-1;
        Aff_tab_2D(T,nbr);
        printf("\tL%d <==> L%d",i+1,etat+1);
        getch();
        /*affichage de la matrice ou cas
        de la réalisation de switch et de retourner la valeur rep
        qu'il va être multiplié par la déterminant*/
    }
    }
    return rep;
}
int Deter_tab_2D(const rational* R,const int nbr,rational pro)
{
    /*Cette fonction est permet de calculer la déterminant*/
    int i;
    for(i=0; i<nbr; i++)
        if(R[i].nmr==0)
        {
            printf("\n\n");
            printf("Determinant= 0 \n(car l'un des elements du diagonale est nulle)");
            return 0;
        }
    printf("\n\n");
    printf("Determinant= %d X ",pro.nmr);
    for(i=0; i<nbr-1; i++)
    {
        pro=mult_r(pro,R[i]);
        if(R[i].dnm==1)
            printf("%d X ",R[i].nmr);
        else printf("%d/%d X ",R[i].nmr,R[i].dnm);
    }
    if(R[nbr-1].dnm==1)
            printf("%d ",R[i].nmr);
    else printf("%d/%d ",R[nbr-1].nmr,R[nbr-1].dnm);
    pro=mult_r(pro,R[nbr-1]);
    pro=redc_r(pro);
    if(pro.dnm==1)
            printf("= %d \n",pro.nmr);
    else printf("= %d/%d \n",pro.nmr,pro.dnm);
    return 0;
}
rational Cramer_det(const rational** T,const int nbr,int k,int j)
{
    /* Cette fonction est une fonction récursive qu'elle prend comme argument
    la matrice "T",l'odre "nbr",l'indice de la linge ignorée "k"
    et l'indice du colonne.
    Et elle retourne double à la fins de la fonction ou bien à
    la vérification de la condition (1)*/
    int i;
    rational det,pos;
    det.dnm=pos.dnm=pos.nmr=1;
    det.nmr=0;
    for(i=0; i<nbr; i++)
    {
        if(nbr==1)
            return T[i][j]; /* (1) return lorsque on a un seul
                            élément dans la matrice*/
        if(!T[i][j].nmr){   // Ignorer les cas de 0
            pos.nmr=pos.nmr*-1;
            continue;
        }
        rational** R;
        R=allo_tab_2D(nbr-1);
        Init_cramer(R,T,i,j,nbr);
        det=add_r(det,mult_r(mult_r(pos,T[i][j]),Cramer_det(R,nbr-1,i,j)));
        pos.nmr=pos.nmr*-1;
        R=deallo_tab_2D(R,nbr-1);
    }
    return det;
    //complixté temporelle: O(2^nbr)
    /*complixté spatiale:
    Sigma des (k*sizeof(rational*)+k*sizeof(rational))
    tels que k= 0->nbr*/

}
void Init_cramer(rational** R,const rational** T,int l,int m,const int nbr)
{
    /* Cette fonction est permet d'initialiser une matrice
     tels que : l'odre de R = (l'odre de T)-1 */
    int i,j,a,b;
    a=b=0;
    for(i=0;i<nbr;i++){
        if(i==l)
            continue;
        for(j=0;j<nbr;j++){
            if(j==m)
                continue;
            R[a][b].nmr=T[i][j].nmr;
            R[a][b].dnm=T[i][j].dnm;
            b++;
        }
    b=0;
    a++;
    }
    //complixié temporelle: O(nbr^2)
    /*complixié spatiale: On a i,j,a,b 4 int donc
    CS=4*sizeof(int)=4*4=16 bytes*/
}
rational** deallo_tab_2D(rational** T,const int nbr)
{
    /*Cette fonction est permet de libèrer la zone mémoire allouee */
    int i;
    for(i=0; i<nbr; i++)
        free(T[i]);
    free(T);
    T=NULL;
    return T;
}
rational* deallo_pivo(rational* R)
{
    /*Cette fonction est permet de libèrer la zone mémoire allouee */
    free(R);
    R=NULL;
    return R;
}
