#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED

#include <stdio.h>
#include <math.h>

//Projection-----------------------------------------------------------------
#define xy                      0
#define xz                      1
#define yz                      2
//-------------------------Constantes-------------------------
#define PI 3.14159
#define MaxNPoints              200     // maximo numero de puntos
#define V_SON                   340     // velocidad del sonido
#define DUR_SIM                 1000    // velocidad simulacion
#define N_REC                   1
using namespace std;

class vector
{
public:
    double x,y,z;
    vector operator+(vector v2){//suma
        vector v1;
        v1.x=x+v2.x;
        v1.y=y+v2.y;
        v1.z=z+v2.z;
        return v1;
    };
    vector operator-(vector v2){//substraccion
        vector v1;
        v1.x=x-v2.x;
        v1.y=y-v2.y;
        v1.z=z-v2.z;
        return v1;
    };
   vector operator*(double f){//multiplicacion por escalar
        vector v;
        v.x=x*f;
        v.y=y*f;
        v.z=z*f;
        return v;
    };
    vector operator/(double f){//division por escalar
        vector v;
        v.x=x/f;
        v.y=y/f;
        v.z=z/f;
        return v;
    };
    double operator*(vector v){//producto escalar
        try{
            double f;
            f=x*v.x+y*v.y+z*v.z;
            return f;
        }
        catch(std::exception& e){
            double f=0;
            return f;
        }
    };
    vector operator/(vector v2){//producto vectorial
        vector v1;
        v1.x=y*v2.z-z*v2.y;
        v1.y=-x*v2.z+z*v2.x;
        v1.z=x*v2.y-y*v2.x;
        return v1;
    };
    void operator*=(double f){//multiplicacion por escalar
        x*=f;
        y*=f;
        z*=f;
    };
    void operator/=(double f){//division por un escalar
        x/=f;
        y/=f;
        z/=f;
    };
    void operator=(double f){//mismo valor para x,y y z
        x=y=z=f;
    };
    bool operator==(vector v){//igualdades entre vectores
        if(x==v.x&&y==v.y&&z==v.z)
            return 1;
        else
            return 0;
    };
    bool operator!=(vector v){//desigualdad entre vectores
        if(x==v.x&&y==v.y&&z==v.z)
            return 0;
        else
            return 1;
    };
    double Angle(vector v){//Angulo entre dos vectores
        double angle,f;
        f=x*v.x+y*v.y+z*v.z;//producto escalar
        f/=sqrt(x*x+y*y+z*z);
        f/=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
        angle =acos(f)/M_PI*180;
        return angle;
    };
    double Max(void){//valor maximo entre las coordenadas
        double max;
        if(x>y){
            max=x;
        } else {
            max=y;
        }
        if(z>max){
            max=z;
        }
        return max;
    };
    double Min(void){//valor minimo entre las coordenadas
        double min;
        if(x<y){
            min=x;
        } else {
            min=y;
        }
        if(z<min){
            min=z;
        }
        return min;
    };
    vector Abs(void)
    {
        //valor absoluto de las coordenadas
        vector v;
        v.x=fabs(x);
        v.y=fabs(y);
        v.z=fabs(z);
        return v;
    }
    double Module(void){
        return sqrt(x*x+y*y+z*z);
    }

    friend std::ostream& operator<<(std::ostream& os, const vector& v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
};
//-------------------------------------------------------------------------------------------------------
class point
{
public:
    double x,y,z;
    point(){
        x=0;
        y=0;
        z=0;
    };
    point(double X, double Y, double Z) {
        x = X;
        y = Y;
        z = Z;
    };
    point operator+(vector v){//traslado
        point p;
        p.x=x+v.x;
        p.y=y+v.y;
        p.z=z+v.z;
        return p;
    };
    point operator+(point p2){//suma dos coordenadas de dos puntos
        point p1;
        p1.x=x+p2.x;
        p1.y=y+p2.y;
        p1.z=z+p2.z;
        return p1;
    };
    vector operator-(point p){//substraccion
        vector v;
        v.x=x-p.x;
        v.y=y-p.y;
        v.z=z-p.z;
        return v;
    };
    point operator*(double f){//multiplicacion por escalar
        point p;
        p.x=x*f;
        p.y=y*f;
        p.z=z*f;
        return p;
    };
    point operator/(double f){//division por escalar
        point p;
        p.x=x/f;
        p.y=y/f;
        p.z=z/f;
        return p;
    };

    void operator=(double f){//mismo valor para x,y y z
        x=y=z=f;
    };
    bool operator==(point p){//igualdades entre puntos
        if(x==p.x&&y==p.y&&z==p.z)
            return 1;
        else
            return 0;
    };
    bool operator!=(point p){//desigualdad entre puntos
        if(x==p.x&&y==p.y&&z==p.z)
            return 0;
        else
            return 1;
    };
    void clear(){
        x=0;
        y=0;
        z=0;
    };
    double Max(void){//valor maximo entre las coordenadas
        double max;
        if(x>y){
            max=x;
        } else {
            max=y;
        }
        if(z>max){
            max=z;
        }
        return max;
    };
    double Min(void){//valor minimo entre las coordenadas
        double min;
        if(x<y){
            min=x;
        } else {
            min=y;
        }
        if(z<min){
            min=z;
        }
        return min;
    };
    point Abs(void){//valor absoluto de las coordenadas
        point v;
        v.x=fabs(x);
        v.y=fabs(y);
        v.z=fabs(z);
        return v;
    };
    double distancia(point p2){
        return sqrt((p2.x-x)*(p2.x-x)+(p2.y-y)*(p2.y-y)+(p2.z-z)*(p2.z-z));
    };
    vector restaPuntos(point b){
        vector v;
        v.x=x-b.x;
        v.y=y-b.y;
        v.z=z-b.z;
        return v;
    };

    friend std::ostream& operator<<(std::ostream& os, const point& p) {
        os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
        return os;
    }
};
//----------------------------------------------------------------------------------
class triangle{
public:
    point p0,p1,p2,bc;  //triangle points
    int Projection;     //projection
    double a0;          //a0 constante para calculos futuros
    int ID;             //identificador unico

    triangle(){
        p0=0;
        p1=0;
        p2=0;
        bc=0;
        Projection=0;
        a0=0;
        ID=0;
    };
    void operator=(triangle t){
        p0=t.p0;
        p1=t.p1;
        p2=t.p2;
        bc=t.bc;
        Projection=t.Projection;
        a0=t.a0;
        ID=t.ID;
    };
    void clear(){
        p0=0;
        p1=0;
        p2=0;
        bc=0;
        Projection=0;
        a0=0;
        ID=0;
    };
    void Centroid(){
        bc=(p0+p1+p2)/3;
    };
    double solidAngle(point b){
        double area=0.0,d=0.2;
        triangle t;
        vector v0,v1,v2;
        v0=p0-b;
        v1=p1-b;
        v2=p2-b;
        v0=v0/v0.Module();
        v1=v1/v1.Module();
        v2=v2/v2.Module();
        t.p0=b+(v0*d);
        t.p1=b+(v1*d);
        t.p2=b+(v2*d);
        area=t.TriangleArea();
        return area;
    };
    double TriangleArea(){
        double a;
        vector v=(p1-p0)/(p2-p0);
        a=0.5*v.Module();
        return a;
    };
    void CalculateProjection() {
        vector n;
        double x0,y0,z0,x1,y1,z1,x2,y2,z2;
        x0=p0.x;
        y0=p0.y;
        z0=p0.z;
        x1=p1.x;
        y1=p1.y;
        z1=p1.z;
        x2=p2.x;
        y2=p2.y;
        z2=p2.z;
        n=(p1-p0)/(p2-p0);
        n.x=n.x*n.x;
        n.y=n.y*n.y;
        n.z=n.z*n.z;
        if((n.x>=n.y)&&(n.x>=n.z)) {                        //projeção yz
            Projection=yz;
            a0=1/(-y1*z0+y2*z0+y0*z1-y2*z1-y0*z2+y1*z2 + 0.000001);
        }
        if((n.y>=n.x)&&(n.y>=n.z)) {                        //projeção xz
            Projection=xz;
            a0=1/(-x1*z0+x2*z0+x0*z1-x2*z1-x0*z2+x1*z2 + 0.000001);
        }
        if((n.z>=n.x)&&(n.z>=n.y)) {                        //projeção xy
            Projection=xy;
            a0=1/(-x1*y0+x2*y0+x0*y1-x2*y1-x0*y2+x1*y2 + 0.000001);
        }
    };

    friend std::ostream& operator<<(std::ostream& os, const triangle& t) {
        os << "Triangle points: ("
           << "p0: (" << t.p0.x << ", " << t.p0.y << ", " << t.p0.z << "), "
           << "p1: (" << t.p1.x << ", " << t.p1.y << ", " << t.p1.z << "), "
           << "p2: (" << t.p2.x << ", " << t.p2.y << ", " << t.p2.z << "))";
        return os;
    }
};
//---------------------------------------------------------------------------
class source
{
public:
    point p;                //Posicion
    vector *Rays;           //Direcciones de partida de los rayos
    int NRAYS;              // Numero de rayos
    double eF;              //Energia inicial de la fuente

    source()
    {
        p=0.0;
        eF=0.0;
        NRAYS=0;
        Rays=NULL;
    };
    void createRays(double NumberOfRays)
    {
        //matriz das Arestas {1o ponto da aresta, 2o ponto da aresta, Posição dos pontos da aresta na matriz Rays}
        int A[30][3]= {{0,1,0}, {0,2,0}, {0,3,0}, {0,4,0}, {0,5,0},
            {1,6,0}, {2,6,0}, {2,7,0}, {3,7,0}, {3,8,0},
            {4,8,0}, {4,9,0}, {5,9,0}, {5,10,0},{1,10,0},
            {6,11,0},{7,11,0},{8,11,0},{9,11,0},{10,11,0},
            {1,2,0}, {2,3,0}, {3,4,0}, {4,5,0}, {5,1,0},
            {6,7,0}, {7,8,0}, {8,9,0}, {9,10,0},{10,6,0}
        };
        //matriz dos triangulos {1a aresta, 2a aresta, [0] V em pé [-1] V de cabeça pra baixo}
        int T[20][3]= {{0,1,0},   {1,2,0},   {2,3,0},   {3,4,0},   {4,0,0},
            {5,6,-1},  {6,7,0},   {7,8,-1},  {8,9,0},   {9,10,-1},
            {10,11,0}, {11,12,-1},{12,13,0}, {13,14,-1},{14,5,0},
            {15,16,-1},{16,17,-1},{17,18,-1},{18,19,-1},{19,15,-1}
        };
        int i,j,k,n,m,RAY;
        double S,R,xB,yB,zB,xC,yC,zC,c[8];
        //create Rays matrix
        //cout<<n;
        if(NRAYS>0)
            delete[] Rays;
        n=int(floor(sqrt((NumberOfRays-2)/10)+0.5));
        NRAYS=int(2+10*pow(n,2));
        Rays=new vector[NRAYS];
        //calculating the icosaedron vertives
        S=2/sqrt(5);
        R=(5-sqrt(5))/5;
        Rays[0].x=0;
        Rays[0].y=0;
        Rays[0].z=1;
        for(i=1; i<6; i++)
        {
            Rays[i].x=S*cos((PI*i*72)/180);
            Rays[i].y=S*sin((PI*i*72)/180);
            Rays[i].z=1-R;
            Rays[i+5].x=S*cos((72*PI*i)/180+(36*PI)/180);
            Rays[i+5].y=S*sin((72*PI*i)/180+(36*PI)/180);
            Rays[i+5].z=R-1;
        }
        Rays[11].x=0;
        Rays[11].y=0;
        Rays[11].z=-1;
        RAY=12;
        //calculating the rays on the icosaedron edges
        for(j=0; j<30; j++)
        {
            A[j][2]=RAY;
            xB=Rays[A[j][0]].x;
            yB=Rays[A[j][0]].y;
            zB=Rays[A[j][0]].z;
            xC=Rays[A[j][1]].x;
            yC=Rays[A[j][1]].y;
            zC=Rays[A[j][1]].z;
            c[0]=pow(xC,2)*(pow(yB,2)+pow(zB,2))+pow(yC*zB-yB*zC,2)-2*xB*xC*(yB*yC+zB*zC)+pow(xB,2)*(pow(yC,2)+pow(zC,2));
            c[1]=acos(xB*xC+yB*yC+zB*zC);
            c[2]=-xC*(yB*yC+zB*zC)+xB*(pow(yC,2)+pow(zC,2));
            c[3]=xC*(pow(yB,2)+pow(zB,2))-xB*(yB*yC+zB*zC);
            c[4] = pow(xC, 2) * yB - xB * xC * yC + zC * (-yC * zB + yB * zC);
            c[5] = -xB * xC * yB + pow(xB, 2) * yC + zB * (yC * zB - yB * zC);
            c[6] = pow(xC, 2) * zB - xB * xC * zC + yC * (yC * zB - yB * zC);
            c[7] = -xB * xC * zB + pow(xB, 2) * zC + yB * (-yC * zB + yB * zC);
            for(i=1; i<n; i++)
            {
                Rays[RAY].x=(c[2]*cos(i*c[1]/n)+c[3]*cos((n-i)*c[1]/n))/c[0];
                Rays[RAY].y=(c[4]*cos(i*c[1]/n)+c[5]*cos((n-i)*c[1]/n))/c[0];
                Rays[RAY].z=(c[6]*cos(i*c[1]/n)+c[7]*cos((n-i)*c[1]/n))/c[0];
                RAY++;
            }
        }
        //calculating the rays on the icosaedron faces
        for(k=0; k<20; k++)
            for(j=1; j<n; j++)
            {
                xB=Rays[A[T[k][0]][2]+j-1].x;
                yB=Rays[A[T[k][0]][2]+j-1].y;
                zB=Rays[A[T[k][0]][2]+j-1].z;
                xC=Rays[A[T[k][1]][2]+j-1].x;
                yC=Rays[A[T[k][1]][2]+j-1].y;
                zC=Rays[A[T[k][1]][2]+j-1].z;
                c[0] = pow(xC,2)*(pow(yB,2)+pow(zB,2))+pow(yC*zB-yB*zC,2)-2*xB*xC*(yB*yC+zB*zC)+pow(xB,2)*(pow(yC,2)+pow(zC,2));
                c[1]=acos(xB*xC+yB*yC+zB*zC);
                c[2]=-xC*(yB*yC+zB*zC)+xB*(pow(yC,2)+pow(zC,2));
                c[3]=xC*(pow(yB,2)+pow(zB,2))-xB*(yB*yC+zB*zC);
                c[4]=pow(xC,2)*yB-xB*xC*yC+zC*(-yC*zB+yB*zC);
                c[5]=-xB*xC*yB+pow(xB,2)*yC+zB*(yC*zB-yB*zC);
                c[6]=pow(xC,2)*zB-xB*xC*zC+yC*(yC*zB-yB*zC);
                c[7]=-xB*xC*zB+pow(xB,2)*zC+yB*(-yC*zB+yB*zC);
                if(T[k][2]==0)m=j;
                else m=n-j;
                for(i=1; i<m; i++)
                {
                    Rays[RAY].x=(c[2]*cos(i*c[1]/m)+c[3]*cos((m-i)*c[1]/m))/c[0];
                    Rays[RAY].y=(c[4]*cos(i*c[1]/m)+c[5]*cos((m-i)*c[1]/m))/c[0];
                    Rays[RAY].z=(c[6]*cos(i*c[1]/m)+c[7]*cos((m-i)*c[1]/m))/c[0];
                    RAY++;
                }
            }
    };
};
//---------------------------------------------------------------------------
class pointGrid {
public:
    point** p;      //Matriz dinámica de puntos
    int I, J;         //Número de puntos

    pointGrid() {
        I = 0;
        J = 0;
        p = NULL;
    };
    void init(int x, int y) {
        I = x;
        J = y;
        p = new point * [I];
        for (int i = 0; i < I; i++) {
            p[i] = new point[J];
            for (int j = 0; j < J; j++)
                p[i][j] = 0.0;
        }
    };
};
//-------------------------------------------------------------------------------------------------------
class plane {
public:
    int         NP;                     //Number of Points
    point* p;                     //plane Points
    int         NT;                     //Number of Triangles
    triangle* t;                     //plane Triangles
    vector      n;                      //Normal Vector

    plane() {
        int P, T;
        NP = 0;
        p = NULL;
        NT = 0;
        t = NULL;
        n = 0;
    }
    void NewTriangle(int N) {
        int T;
        triangle* tt;
        tt = new triangle[NT + N];
        for (T = 0; T < NT; T++) {
            tt[T] = t[T];
        }
        for (T = NT; T < NT + N; T++) {
            tt[T].clear();
        }
        if (NP > 0) {
            delete[] t;
            t = NULL;
        }
        t = tt;
        NT += N;
    };
    void NewPoints(int N) {
        int P;
        point* tp;
        tp = new point[NP + N];
        for (P = 0; P < NP; P++) {
            tp[P] = p[P];
        }
        for (P = NP; P < NP + N; P++) {
            tp[P].clear();
        }
        if (NP > 0) {
            delete[] p;
            p = NULL;
        }
        p = tp;
        NP += N;
    };
    void PointGenTriangle() { //Genera 2 triangulos a partir de los vertices de un cuadrado
        NewTriangle(NP - 2);
        int i = 1;
        for (int T = 0; T < NT; T++) {
            i--;
            t[T].p0.x = p[i].x;
            t[T].p0.y = p[i].y;
            t[T].p0.z = p[i].z;
            i++;
            if (i == NP) i = 0;
            t[T].p1.x = p[i].x;
            t[T].p1.y = p[i].y;
            t[T].p1.z = p[i].z;
            i++;
            if (i == NP) i = 0;
            t[T].p2.x = p[i].x;
            t[T].p2.y = p[i].y;
            t[T].p2.z = p[i].z;
            i++;
        }
    };
    void MoreTriangle(int nd) { //Genera m�s tri�ngulos a partir de una malla con nd divisiones
        if (NP == 4) {
            int i, j, cont;   //Contadores
            pointGrid np;   //Matriz din�mica de puntos
            vector v1, v2;  //Vectores directores en cada lado del cuadrado
            double m1, m2;  //m�dulos de los Vectores directores
            double p1, p2;  //tama�o del paso
            v1 = p[1] - p[0];   //Vector director 1 dgenerado por el v�rtice 1 y 2
            v2 = p[2] - p[1];   //Vector director 2 generado por el v�rtice 2 y 3
            m1 = v1.Module(); //m�dulo del Vector director 1
            m2 = v2.Module(); //m�dulo del Vector director 2
            v1 = v1 / m1;       //Vector director 1 unitario
            v2 = v2 / m2;       //Vector director 2 unitario
            p1 = m1 / nd;       //paso 1
            p2 = m2 / nd;       //paso 2

            np.init(nd + 1, nd + 1);//Lleno la matriz de puntos, seg�n los v�rtices del cuadrado inicial.
            for (i = 0; i <= nd; i++) {
                np.p[i][0] = p[0] + (v1 * (p1 * i));
                for (j = 1; j <= nd; j++)
                    np.p[i][j] = np.p[i][0] + (v2 * (p2 * j));
            }

            plane* a_p = new plane[nd * nd];
            cont = 0;
            for (i = 0; i < nd; i++) {
                for (j = 0; j < nd; j++) {
                    a_p[cont].clear();
                    a_p[cont].NewPoints(4);
                    a_p[cont].p[0] = np.p[i][j];
                    a_p[cont].p[1] = np.p[i + 1][j];
                    a_p[cont].p[2] = np.p[i + 1][j + 1];
                    a_p[cont].p[3] = np.p[i][j + 1];
                    a_p[cont].PointGenTriangle();
                    cont++;
                }
            }

            cont = 0;
            NewTriangle(2 * nd * nd);
            for (int i = 0; i < nd * nd; i++) {
                for (int j = 0; j < a_p[i].NT; j++) {
                    t[cont] = a_p[i].t[j];
                    cont++;
                }
            }
            delete a_p;
            a_p = NULL;
        }
    };
    void clear() {
        NP = 0;
        delete[] p;
        p = NULL;
        NT = 0;
        delete[] t;
        t = NULL;
        n = 0;
    };
    /*friend std::ostream& operator<<(std::ostream& os, const plane& pl) {
        os << "Plane vertices: (" << pl.p[0] << ", " << pl.p[1] << ", " << pl.p[2] << ")";
        return os;
    }*/
    void PrintPlane(const plane& pl)
    {
        cout << "Plane vertices: (" << pl.p[0] << ", " << pl.p[1] << ", " << pl.p[2] << ")" << endl;
    }
};
//---------------------------------------------------------------------------
class receptor {
public:
    point p;                //Posición
    double radioReceptor;   //Radio del receptor
    double* eR;             //Energía recibida en el receptor
    int NIt;                //Instantes de tiempo considerados
    int Nrays;             // numero de rayos que chocaron

    receptor() {
        p = 0.0;
        eR = NULL;
        NIt = 0;
        radioReceptor = 0.3;
    };

    void clear() {
        if (NIt > 0) {
            delete eR;
            eR = NULL;
            NIt = 0;
        }
    };

    void createTimeSamples(int durSim) {
        clear();
        NIt = durSim;
        eR = new double[NIt];
        for (int i = 0; i < NIt; i++) {
            eR[i] = 0.0;
        }
    };

    double IntersectionDistance(vector n, point p, vector u, point o) {
        double d, m;
        m = n * u;
        if (m == 0)
            return -1;
        d = (n * (p - o)) / m;
        return d;
    };

    double Module(vector v) {
        double m;
        m = sqrt(v * v);
        return m;
    };

    vector Normal(vector v1) {
        double m;
        vector v2;
        m = Module(v1);
        if (m == 0)
            v2 = 0;
        else
            v2 = v1 / m;
        return v2;

    };

    double solidAngle(point b) {
        double area = 0.0, d = 0.0, r = 0.0, d2 = 0.2;
        d = Module(p - b);
        r=(radioReceptor*d2)/d;
        area = PI * r * r;
        return area;
    };

    void CalculoIncidencia(vector in, point o, double Ener, int t,double maxd)
    {
        // Distancia centro baricentro
        double distanciaEner;
        int temp;
        // Vector normal invertido
        vector normal,u;
        u=Normal(in);
        normal = (Normal(u) * (-1));


        // Calcular punto de incidencia
        double distancia = IntersectionDistance(normal,p,u,o);


        if (distancia > 0 ) {//&& distancia < maxd
            // Punto de incidencia
            point pi = o + (in * distancia);

            // Distancia para comparar
            double discomp = p.distancia(pi);

            // Si es que la distancia
            temp = 0;
            if (discomp <= radioReceptor)
            {
                // Calcular el tiempo del rayo en chocar con el receptor
                distanciaEner = o.distancia(p);
                temp = int((distanciaEner / (V_SON)) * DUR_SIM);
                temp += t;
                // Registrar la energia
                if (temp >= 1000){
                    temp = 999;
                }
                eR[temp] += Ener;
                // Aumentar numero de rayos que4 chocaron
                Nrays++;
            }
        }
    }
};
//---------------------------------------------------------------------------
struct reflection {                     // respuesta al proceso de trazado de rayos.
    point r[MaxNPoints];                // puntos de aplicacion
    double d[MaxNPoints];               // distancia entre punto y punto
    int idTriangle[MaxNPoints];         // id unico del triangulo por cuarto donde se choco
    int Plane[MaxNPoints];              // nro. del plano por cuarto donde se choco
    int Triangle[MaxNPoints];           // nro. del triangulo por plano donde se choco
    int N;                              // numero de reflexiones
    bool lost;                          // Si lost es igual a true, es un reflejo de un rayo perdido.
    int Ray;                            // Número de RAY en vista previa
};
//---------------------------------------------------------------------------
class room {
public:
    int NP;		        // Numero de planos
    plane* p;		    // Planos
    double maxd;	    // Maxima distancia entre dos puntos en la sala.
    int NR;             // Number of receptors
    receptor* rec;        //Number of receptors (microphones)

    room() {
        NP = 0;
        p = NULL;
        NR = 0;
        maxd = 0.0;
        rec=NULL;
    };

    void clear() {
        if (NP > 0) {
            for (int i = 0; i < NP; i++) {
                p[i].clear();
            }
            delete[] p;
            p = NULL;
        }
        NP = 0;
        if (NR > 0) {
            for (int i = 0; i < NR; i++) {
                rec[i].clear();
            }
            delete[] rec;
            rec = NULL;
        }
        NR = 0;
        maxd = 0.0;
    };
    void Init() {
        NP = 0;
        p = NULL;
        NR = 0;
        maxd = 0.0;
        rec=NULL;
    };
    bool Inner(point p, triangle t) { // verificar si un punto pertenece al triangulo
        double a1, a2, x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2;

        x = p.x;
        y = p.y;
        z = p.z;

        x0 = t.p0.x;
        y0 = t.p0.y;
        z0 = t.p0.z;
        x1 = t.p1.x;
        y1 = t.p1.y;
        z1 = t.p1.z;
        x2 = t.p2.x;
        y2 = t.p2.y;
        z2 = t.p2.z;

        if (t.Projection == yz) {                              // proyeccion yz
            a1 = -t.a0 * (-y0 * z + y2 * z + y * z0 - y2 * z0 - y * z2 + y0 * z2);
            a2 = -t.a0 * (y0 * z - y1 * z - y * z0 + y1 * z0 + y * z1 - y0 * z1);
        }
        if (t.Projection == xz) {                              // proyeccion xz
            a1 = -t.a0 * (-x0 * z + x2 * z + x * z0 - x2 * z0 - x * z2 + x0 * z2);
            a2 = -t.a0 * (x0 * z - x1 * z - x * z0 + x1 * z0 + x * z1 - x0 * z1);
        }
        if (t.Projection == xy) {                              // proyeccion xy
            a1 = -t.a0 * (-x2 * y0 + x * y0 + x0 * y2 - x * y2 - x0 * y + x2 * y);
            a2 = t.a0 * (-x1 * y0 + x * y0 + x0 * y1 - x * y1 - x0 * y + x1 * y);
        }

        if ((a1 + a2 <= 1) && (a1 >= 0) && (a2 >= 0))
            return true;
        else
            return false;

    };
    double IntersectionDistance(vector n, point p, vector u, point o) {
        double d, m;
        m = n * u;
        if (m == 0)
            return -1;
        d = (n * (p - o)) / m;
        return d;
    };
    void MaxDistance() {
        maxd = 0;
        float tmpd = 0;
        for (int i1 = 0; i1 < NP; i1++) {
            for (int j1 = 0; j1 < p[i1].NP; j1++) {
                for (int i2 = 0; i2 < NP; i2++) {
                    for (int j2 = 0; j2 < p[i2].NP; j2++) {
                        tmpd = p[i1].p[j1].distancia(p[i2].p[j2]);
                        if (maxd < tmpd)
                            maxd = tmpd;
                    }
                }
            }
        }
    };
    void NewPlanes(int N) {
        int P;
        plane* tp;
        tp = new plane[NP + N];
        for (P = 0; P < NP; P++) {
            tp[P] = p[P];
        }
        for (P = NP; P < NP + N; P++) {
            tp[P].clear();
        }
        if (NP > 0) {
            delete[] p;
            p = NULL;
        }
        p = tp;
        NP += N;
    };
    void NewReceptor(int N) {
        int R;
        receptor* tr;
        tr = new receptor[NR + N];
        for (R = 0; R < NR; R++) {
            tr[R] = rec[R];
        }
        for (R = NR; R < NR + N; R++) {
            tr[R].clear();
        }
        if (NR > 0) {
            delete[] rec;
            rec = NULL;
        }
        rec = tr;
        NR += N;
    };
    double Module(vector v) {
        double m;
        m = sqrt(v * v);
        return m;
    };
    vector Normal(vector v1) {

        double m;
        vector v2;
        m = Module(v1);
        if (m == 0)
            v2 = 0;
        else
            v2 = v1 / m;
        return v2;

    };
    reflection* RayTracing(point origen, vector* Rays, int NRAYS) {

        reflection* ry;
        ry = NULL;

        int
            IntersectedPlane,       // indice del plano interesectado
            IntersectedTriangle,    // indice del triangulo intersectado por plano
            IntersectedTriangleId,  // id unico del tri�ngulo intersectado
            NReflections,           // numero actual de reflexion
            TNReflections,          // numero total de reflexiones
            LostRays = 0;           // contador de rayos perdidos

        double      // distancia al punto de intersecci�n
            d1,
            d2;

        point       // puntos para determinar donde existe intersecci�n con el plano
            p1,
            p2,
            p3;

        bool Stop;  // bandera para detener el procedimiento
        vector v;   // vector incidente
        point o;    // punto de partida (origen del Vector incidente)

        ry = new reflection[NRAYS];

        NReflections = 0;
        TNReflections = 0;

        // inicio del RAY-TRACING
        MaxDistance();
        for (int R = 0; R < NRAYS; R++) {
            v = Rays[R];  // asigno el primer rayo del arreglo original a v
            //cout<<endl<<"Rayo original "<<R+1<<endl;
           //cout<<"("<<v.x<<","<<v.y<<","<<v.z<<")"<<endl;
            o = origen;

            // como no existe aun ninguna reflexion, inicializo con el valor -1
            IntersectedPlane = -1;
            IntersectedTriangle = -1;
            IntersectedTriangleId = -1;

            Stop = false;

            ry[R].N = 0;            // numero de reflexion, inincialmente 0
            ry[R].r[0] = o;         // punto de partida, inicialmente el centro de la fuente
            ry[R].d[0] = 0.0;       // distancia recorrida, inicialmente 0
            ry[R].lost = false;     // reflexion perdida? inicialmente falso
            ry[R].Ray = R;          // rayo asociado a la reflexion

            //cout<<"Reflexiones del rayo "<<R+1<<endl;
            while (!Stop) {   // lazo para realizar varias reflexiones
                d1 = maxd;    // asigno a d1 la m�xima distancia entre puntos de la sala
                for (int P = 0; P < NP; P++) { // recorro todos los planos de la sala
                    if ((p[P].n * v) < 0) { // existe ruta de intersecci�n entre un plano y un Vector //  && (P!=LastIntersectedPlane)
                        d2 = IntersectionDistance(p[P].n, p[P].p[0], v, o);
                        if ((d2 > 0.0) && (d2 < d1)) { // verifico que la distancia de vuelo del rayo no sea cero && que no sea mayor a d1 (m�xima distancia de vuelo o otra ruta de intersecci�n con otro plano)
                            p2 = o + (v * d2); // obtengo el punto de incidencia en el plano
                            for (int T = 0; T < p[P].NT; T++) { // recorro todos los tri�ngulos del plano
                                if (Inner(p2, p[P].t[T])) {// verifica si el punto pertenece al tri�ngulo
                                    // registro la distancia y el punto de intersecci�n
                                    d1 = d2;
                                    p1 = p2;
                                    // registro los identificadores del elemento interesectado
                                    IntersectedPlane = P;
                                    IntersectedTriangle = T;
                                    IntersectedTriangleId = p[P].t[T].ID;
                                    T = p[P].NT; // para forzar la finalizaci�n del recorrido de tri�ngulos
                                }
                            }
                        }
                    }
                }
                if (d1 < maxd && IntersectedPlane != -1) {// si hubo intersecci�n
                    // calculo el Vector reflejo
                    p3 = o;
                    o = p1;// nuevo punto de partida del Vector reflejado
                    v = Normal(v - (p[IntersectedPlane].n * (v * p[IntersectedPlane].n * 2))); // formula del Vector reflejo
                    //cout<<"  Reflexion "<<NReflections<<":"<<"("<<v.x<<","<<v.y<<","<<v.z<<")"<<endl;
                    NReflections++;
                    TNReflections += NReflections;
                    ry[R].r[NReflections] = p1;                                   // puntos de aplicaci�n
                    ry[R].d[NReflections] = Module(p1 - p3);                       // distancia entre punto y punto
                    ry[R].idTriangle[NReflections] = IntersectedTriangleId;       // id unico del tri�ngulo por cuarto donde se choc�
                    ry[R].Plane[NReflections] = IntersectedPlane;                 // nro. del plano por cuarto donde se choc�
                    ry[R].Triangle[NReflections] = IntersectedTriangle;           // nro. del tri�ngulo por plano donde se choc�
                    ry[R].N = NReflections;                                       // Number of reflections

                    // maximo 50 reflexiones
                    if (NReflections > 50) {
                        Stop = true;                  // no realizo mas reflexiones con este rayo
                        NReflections = 0;             // reseteo el contador de reflexiones para el siguiente rayo
                        IntersectedPlane = -1;        // reseteo el identificador de plano intersecado
                        IntersectedTriangle = -1;     // reseteo el identificador de triangulo intersecado
                        IntersectedTriangleId = -1;   // reseteo el identificador �nico de triangulo intersecado
                    }
                }
                else {// no hubo intersecci�n
                    NReflections++;
                    ry[R].lost = true;// reflexion perdida
                    p3 = o + (v * maxd); // define un punto fuera de la sala para ilustrar el rayo perdido
                    ry[R].r[NReflections] = p3;// defino un punto fuera de la sala para ilustrar el rayo perdido
                    ry[R].N = NReflections;// punto fuera de la sala
                    LostRays++;                 // contador de rayos perdidos
                    Stop = true;                  // no realizo mas reflexiones con el rayo perdido
                    NReflections = 0;             // reseteo el contador de reflexiones para el siguiente rayo
                    IntersectedPlane = -1;        // reseteo el identificador de plano intersecado
                    IntersectedTriangle = -1;     // reseteo el identificador de triangulo intersecado
                    IntersectedTriangleId = -1;   // reseteo el identificador �nico de triangulo intersecado
                }
            }
        }
        return ry;
    };
};
//---------------------------------------------------------------------------
class matInt {
public:
    int** m;      //Matriz din�mica de puntos
    int I, J;         //N�mero de puntos

    matInt() {
        I = 0;
        J = 0;
        m = NULL;
    };
    void clear(){
        I = 0;
        J = 0;
        m = NULL;
        delete[] m;
    }
    void init(int x, int y) {
        I = x;
        J = y;
        m = new int* [I];
        for (int i = 0; i < I; i++) {
            m[i] = new int[J];
            for (int j = 0; j < J; j++)
                m[i][j] = 0.0;
        }
    };
};
//-------------------------------------------------------------------------------------------------------
class matDouble {
public:
    double** d;      //Matriz din�mica de doubles
    int I, J;         //N�mero de filas y columnas

    matDouble() {
        I = 0;
        J = 0;
        d = NULL;
    };
    void init(int x, int y) {
        I = x;
        J = y;
        d = new double* [I];
        for (int i = 0; i < I; i++) {
            d[i] = new double[J];
            for (int j = 0; j < J; j++)
                d[i][j] = 0.0;
        }
    };
    void clear() {
        I = 0;
        J = 0;
        delete[] d;
        d = NULL;
    };
};

inline vector VectorUnitario(vector v1) {
    //Calculo de la normal acorde al vector obtenido del producto vectorial
    vector v2;
    double mod = sqrt(v1 * v1);
    if (mod == 0)
        v2 = 0;
    else
        v2 = v1 / mod;
    return v2;
};
inline vector TriangleNormal(triangle t) {
    vector n;
    n = VectorUnitario((t.p1 - t.p0) / (t.p2 - t.p0));//Calculo del producto vectorial de los puntos
    return n;
};

#endif // CLASSES_H_INCLUDED
