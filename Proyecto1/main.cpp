#include <iostream>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "classes.h"

// Definiciones
#define PI 3.14159

// Variables globales
int numTri = 0;
room r;
source s;
point origin = point(1.5f, 1.5f, 1.5f);//Coordenadas de la fuente
reflection* reflexiones;
int nDivs = 2;
;
/*MATRICES PARA EL ROOM*/
matInt matRoomTime; //Matriz de tiempo de Triangulos de la sala
matDouble matRoomDist; //Matriz de distancia de Triangulos de la sala
matDouble matRoomAngles; // Matriz de angulos de Triangulos de la sala
matDouble MET; // Matriz de energía

/*MATRICES PARA EL RECEPTOR*/
matInt matRecTime; //Matriz de tiempo de Triangulos del receptor
matDouble matRecDist; //Matriz de distancia de Triangulos del receptor
matDouble matRecAngles; // Matriz de angulos de Triangulos del receptor

// Funciones
void generarCubo();
void calcular_eneRes();
void guardar_archivo(const std::string& nombre, std::ofstream& archivoCSV);

using namespace std;

int main()
{
    generarCubo();
    calcular_eneRes();

    //Guardar la matriz de espacio temporal (MET) en un archivos csv
    cout<<"********Guardando matriz MET en el archivo********"<<endl;
    ofstream archivoCSV;
    guardar_archivo("datos.csv", archivoCSV);
    for (int i = 0; i < numTri; ++i) {
        for (int j = 0; j < DUR_SIM; ++j) {
            archivoCSV << MET.d[i][j] << " ";
        }
        archivoCSV << endl;
    }
    archivoCSV.close();

    return 0;
}

void guardar_archivo(const string& nombre, ofstream& archivoCSV) {
    archivoCSV.open(nombre);

    if (!archivoCSV.is_open()) {
        std::cerr << "Error al abrir el archivo CSV." << std::endl;
    }
}

void generarCubo(){
    //Inicializaciones cubo
    r.clear();
    //Inicializaciones cubo
    r.NewPlanes(6);
    //-------------square 1 back
    r.p[0].NewPoints(4);
    r.p[0].p[0].x=-2.0f;
    r.p[0].p[0].y=2.0f;
    r.p[0].p[0].z=2.0f;
    r.p[0].p[1].x=-2.0f;
    r.p[0].p[1].y=-2.0f;
    r.p[0].p[1].z=2.0f;
    r.p[0].p[2].x=-2.0f;
    r.p[0].p[2].y=-2.0f;
    r.p[0].p[2].z=-2.0f;
    r.p[0].p[3].x=-2.0f;
    r.p[0].p[3].y=2.0f;
    r.p[0].p[3].z=-2.0f;
    r.p[0].MoreTriangle(nDivs);
    //r.p[0].PointGenTriangle();

    //-------------square 2 front
    r.p[1].NewPoints(4);
    r.p[1].p[0].x=2.0f;
    r.p[1].p[0].y=2.0f;
    r.p[1].p[0].z=2.0f;
    r.p[1].p[1].x=2.0f;
    r.p[1].p[1].y=2.0f;
    r.p[1].p[1].z=-2.0f;
    r.p[1].p[2].x=2.0f;
    r.p[1].p[2].y=-2.0f;
    r.p[1].p[2].z=-2.0f;
    r.p[1].p[3].x=2.0f;
    r.p[1].p[3].y=-2.0f;
    r.p[1].p[3].z=2.0f;
    r.p[1].MoreTriangle(nDivs);
    //r.p[1].PointGenTriangle();

    //-------------square 3 left
    r.p[2].NewPoints(4);
    r.p[2].p[0].x=-2.0f;
    r.p[2].p[0].y=-2.0f;
    r.p[2].p[0].z=2.0f;
    r.p[2].p[1].x=2.0f;
    r.p[2].p[1].y=-2.0f;
    r.p[2].p[1].z=2.0f;
    r.p[2].p[2].x=2.0f;
    r.p[2].p[2].y=-2.0f;
    r.p[2].p[2].z=-2.0f;
    r.p[2].p[3].x=-2.0f;
    r.p[2].p[3].y=-2.0f;
    r.p[2].p[3].z=-2.0f;
    r.p[2].MoreTriangle(nDivs);
    //r.p[2].PointGenTriangle();

    //-------------square 4 right
    r.p[3].NewPoints(4);
    r.p[3].p[0].x=2.0f;
    r.p[3].p[0].y=2.0f;
    r.p[3].p[0].z=2.0f;
    r.p[3].p[1].x=-2.0f;
    r.p[3].p[1].y=2.0f;
    r.p[3].p[1].z=2.0f;
    r.p[3].p[2].x=-2.0f;
    r.p[3].p[2].y=2.0f;
    r.p[3].p[2].z=-2.0f;
    r.p[3].p[3].x=2.0f;
    r.p[3].p[3].y=2.0f;
    r.p[3].p[3].z=-2.0f;
    r.p[3].MoreTriangle(nDivs);
    //r.p[3].PointGenTriangle();

    //-------------square 5 top
    r.p[4].NewPoints(4);
    r.p[4].p[0].x=-2.0f;
    r.p[4].p[0].y=-2.0f;
    r.p[4].p[0].z=2.0f;
    r.p[4].p[1].x=-2.0f;
    r.p[4].p[1].y=2.0f;
    r.p[4].p[1].z=2.0f;
    r.p[4].p[2].x=2.0f;
    r.p[4].p[2].y=2.0f;
    r.p[4].p[2].z=2.0f;
    r.p[4].p[3].x=2.0f;
    r.p[4].p[3].y=-2.0f;
    r.p[4].p[3].z=2.0f;
    r.p[4].MoreTriangle(nDivs);
    //r.p[4].PointGenTriangle();

    //-------------square 1 bottom
    r.p[5].NewPoints(4);
    r.p[5].p[0].x=-2.0f;
    r.p[5].p[0].y=2.0f;
    r.p[5].p[0].z=-2.0f;
    r.p[5].p[1].x=-2.0f;
    r.p[5].p[1].y=-2.0f;
    r.p[5].p[1].z=-2.0f;
    r.p[5].p[2].x=2.0f;
    r.p[5].p[2].y=-2.0f;
    r.p[5].p[2].z=-2.0f;
    r.p[5].p[3].x=2.0f;
    r.p[5].p[3].y=2.0f;
    r.p[5].p[3].z=-2.0f;
    r.p[5].MoreTriangle(nDivs);
    //r.p[5].PointGenTriangle();

    for (int i = 0; i < r.NP; i++)
    {
        cout << "-------------------------------------------------------\n";
        cout << "Vértices del plano " << i+1 << ":\n";
        for (int j = 0; j < r.p[i].NP; j++)
        {
            cout << "  Vértice " << j << ": (" << r.p[i].p[j].x << ", " << r.p[i].p[j].y << ", " << r.p[i].p[j].z << ")\n";
        }
        cout << "Triangulos generados del plano: " << i << ":\n";
        for (int k = 0; k < r.p[i].NT; k++)
        {
            cout << "  Triángulo " << k << ":\n";
            cout << "    Punto 0: (" << r.p[i].t[k].p0.x << ", " << r.p[i].t[k].p0.y << ", " << r.p[i].t[k].p0.z << ")\n";
            cout << "    Punto 1: (" << r.p[i].t[k].p1.x << ", " << r.p[i].t[k].p1.y << ", " << r.p[i].t[k].p1.z << ")\n";
            cout << "    Punto 2: (" << r.p[i].t[k].p2.x << ", " << r.p[i].t[k].p2.y << ", " << r.p[i].t[k].p2.z << ")\n";
        }
    }

    //double areaTriangulo=r.p[0].t[0].TriangleArea();
    //cout<<"Area triangulo "<<areaTriangulo<<endl;

    cout << "-------------------------------------------------------\n";
    cout <<endl<< "*******************************************************\n";
    //Inicializar normales de los planos
    // se calculan las normales con la normal de su primer triangulo
    // se calcula el baricentro de los triángulos
    int cont_t=0;
    cout<<"Baricentro de cada triangulo:\n";
    // Creación de baricentros de triángulos de los planos
    for (int i=0; i<r.NP; i++)
    {
        r.p[i].n=TriangleNormal(r.p[i].t[0]);//Calculo de la normal del plano
        cout<<endl<<"Plano "<<i+1<<":"<<endl;
        for (int j=0; j<r.p[i].NT; j++)
        {
            r.p[i].t[j].CalculateProjection();
            r.p[i].t[j].Centroid();
            r.p[i].t[j].ID=cont_t;
            cont_t++;
            cout<<"Baricentro triangulo "<<j<<": "<<r.p[i].t[j].bc<<endl;
        }
    }
    // CREACIÓN DE RECEPTORES
    r.NewReceptor(N_REC);
    r.rec[0].p.x=-1.5;
    r.rec[0].p.y=-1.5;
    r.rec[0].p.z=-1.5;
    cout<<endl<<"Posiciones de cada receptor:\n";
    for(int i=0;i<N_REC;i++){
        cout<<"Posicion receptor "<<i<<": "<<r.rec[i].p<<endl;
    }

    // Creación de matrices de tiempo, distancia y porcentajes para el room
    numTri=cont_t;
    matRoomTime.init(numTri,numTri);
    matRoomDist.init(numTri,numTri);
    matRoomAngles.init(numTri,numTri);

    // Creación de matrices de tiempo, distancia y porcentajes para los receptores
    matRecTime.init(N_REC, numTri);
    matRecDist.init(N_REC, numTri);
    matRecAngles.init(N_REC, numTri);

    //Arreglo suma de area total ( esto servira obtener el porcentaje de energia)
    double* areaT;
    areaT = NULL;
    areaT = new double[numTri];
    int idTri1 = 0;
    int idTri2 = 0;

    for (int i = 0; i < numTri; i++)
    {
        areaT[i] = 0.0;
    }

    //Calculo de distancias, tiempo de vuelo y porcentajes
    for (int i = 0; i < r.NP; i++)
    {
        for (int j = 0; j < r.p[i].NT; j++)
        {
            // Triangulo 1
            idTri1 = r.p[i].t[j].ID;
            double sumafila=0.0;
            for (int k = 0; k < r.NP; k++)
            {
                for (int l = 0; l < r.p[k].NT; l++)
                {
                    // Triangulo 2
                    idTri2 = r.p[k].t[l].ID;
                    if (i != k)
                    {
                        matRoomDist.d[idTri1][idTri2] = r.p[i].t[j].bc.distancia(r.p[k].t[l].bc);
                        matRoomTime.m[idTri1][idTri2] = int(1000 * matRoomDist.d[idTri1][idTri2] / V_SON);
                        matRoomAngles.d[idTri1][idTri2] = r.p[k].t[l].solidAngle(r.p[i].t[j].bc);
                        sumafila+= matRoomAngles.d[idTri1][idTri2];
                    }
                }
            }
            areaT[idTri1]=sumafila;

            // Cálculo de las distancias, tiempo de vuelo y porcentajes de receptores
            for (int m = 0; m < r.NR; m++) {
                matRecDist.d[m][idTri1] = r.rec[m].p.distancia(r.p[i].t[j].bc);
                matRecTime.m[m][idTri1] = int(1000 * matRecDist.d[m][idTri1] / V_SON);
                matRecAngles.d[m][idTri1] = r.rec[m].solidAngle(r.p[i].t[j].bc);
            }
        }
    }

    // MATRIZ DE PORCENTAJES
    //Sala
    for (int i = 0; i < numTri; i++) {
        for (int j = 0; j< numTri; j++) {
            matRoomAngles.d[i][j] = matRoomAngles.d[i][j]/areaT[i]; //Normalización de de ángulos solidos de la Sala
        };
    }
    //Receptores
    for(int i=0;i<r.NR;i++){
        /*double sumaFilaR=0.0;
        for(int j=0;j<numTri;j++){
            sumaFilaR+=matRecAngles.d[i][j];
        }*/
        for(int k=0;k<numTri;k++){
            //matRecAngles.d[i][k]=matRecAngles.d[i][k]/areaT[i];
            matRecAngles.d[i][k]=matRecAngles.d[i][k]/areaT[k];
        }
    }
    cout<<"********Guardando matriz matRoomDist en el archivo********"<<endl;
    ofstream archivoCSV1;
    ofstream archivoCSV2;
    ofstream archivoCSV3;
    guardar_archivo("matRoomDist.csv", archivoCSV1);
    for (int i = 0; i < numTri; ++i) {
        for (int j = 0; j < numTri; ++j) {
            archivoCSV1 << matRoomDist.d[i][j] << " ";
        }
        archivoCSV1 << endl;
    }
    cout<<"********Guardando matriz matRoomTime en el archivo********"<<endl;
    guardar_archivo("matRoomTime.csv", archivoCSV2);
    for (int i = 0; i < numTri; ++i) {
        for (int j = 0; j < numTri; ++j) {
            archivoCSV2 << matRoomTime.m[i][j] << " ";
        }
        archivoCSV2 << endl;
    }
    cout<<"********Guardando matriz matRoomAngles en el archivo********"<<endl;
    guardar_archivo("matRoomAngles.csv", archivoCSV3);
    for (int i = 0; i < numTri; ++i) {
        for (int j = 0; j < numTri; ++j) {
            archivoCSV3 << matRoomAngles.d[i][j] << " ";
        }
        archivoCSV3 << endl;
    }

    cout<<"********Guardando matriz matRecDist en el archivo********"<<endl;
    ofstream archivoCSV4;
    ofstream archivoCSV5;
    ofstream archivoCSV6;
    guardar_archivo("matRecDist.csv", archivoCSV4);
    for (int i = 0; i < r.NR; ++i) {
        for (int j = 0; j < numTri; ++j) {
            archivoCSV4 << matRecDist.d[i][j] << " ";
        }
        archivoCSV4 << endl;
    }
    cout<<"********Guardando matriz matRecTime en el archivo********"<<endl;
    guardar_archivo("matRecTime.csv", archivoCSV5);
    for (int i = 0; i < r.NR; ++i) {
        for (int j = 0; j < numTri; ++j) {
            archivoCSV5 << matRecTime.m[i][j] << " ";
        }
        archivoCSV5 << endl;
    }
    cout<<"********Guardando matriz matRecAngles en el archivo********"<<endl;
    guardar_archivo("matRecAngles.csv", archivoCSV6);
    for (int i = 0; i < r.NR; ++i) {
        for (int j = 0; j < numTri; ++j) {
            archivoCSV6 << matRecAngles.d[i][j] << " ";
        }
        archivoCSV6 << endl;
    }
}


//Calcular energias
void calcular_eneRes(){
    int numRay=12;
    s.createRays(numRay);
    s.eF=120; //Energia de la fuente

    cout<<endl<<"Los rayos generados son: "<<endl;
    for(int i=0;i<s.NRAYS;i++){
        cout<<"Rayo "<<i+1<<": ("<<s.Rays[i].x<<","<<s.Rays[i].y<<","<<s.Rays[i].z<<")"<<endl;
    }

    double eneRay, eneRes; //Energía del rayo y energía residual;
    eneRay = s.eF / s.NRAYS; // Energía Inicial

    //Coeficientes de absorción y difusión
    double alfa, delta;
    alfa = 0.2;
    delta = 0.2;

    reflexiones = NULL; // reflexion
    reflexiones = r.RayTracing(origin, s.Rays, s.NRAYS); // trazado de rayos

    //creacion matriz receptores
    for (int i = 0; i < r.NR; i++)
    {
        r.rec[i].createTimeSamples(DUR_SIM);
    }

    //Asignar la energia en el receptor

    point partida;
    vector vin;
    double t;
    int time;
    vector inicial;
    for (int i = 0; i < r.NR; i++) // Receptores
    {
        for (int j = 0; j < numRay; j++) // Rayos
        {
            // Punto de partida del rayo
            partida.x = s.Rays[j].x;
            partida.y = s.Rays[j].y;
            partida.z = s.Rays[j].z;

            eneRes = eneRay; // Energia inicial del rayo
            t = 0.0;         // Inicio el tiempo en cero para cada rayo

            for (int k = 0; k < reflexiones[j].N; k++) // Numero de reflexiones
            {
                // Actualizacion de energia
                if (k == 0)
                {
                    vin = partida - reflexiones[j].r[k];
                    t += (reflexiones[j].d[k]) / V_SON;
                }
                else
                {
                    vin = reflexiones[j].r[k - 1] - reflexiones[j].r[k];
                    t += (reflexiones[j].d[k]) / V_SON;
                }
                // Busqueda de un choque.
                time = int((t)*DUR_SIM);
                r.rec[i].CalculoIncidencia(vin, reflexiones[j].r[k], eneRes, time,r.maxd);
                eneRes = eneRes * (1 - alfa) * (1 - delta);
            }
        }
    }

    // Creación de matriz energia Room
    MET.init(numTri,DUR_SIM);//asignar espacio en memoria para la matriz espacio temporal

    // Primera asignacion de energia
    double distancia;
    for(int j=0; j<s.NRAYS; j++)
    {
        distancia = 0.0;
        for(int i=0; i<reflexiones[j].N; i++)
        {
            int tri,tim;
            tri=reflexiones[j].idTriangle[i];
            distancia += reflexiones[j].d[i];
            tim=int((distancia/V_SON)*1000);

            if(tim<DUR_SIM)
            {
                //Asignar la energia en matriz espacio temporal
                MET.d[tri][tim]=eneRay*pow(((1-alfa)*delta),i);
            }
        }
    }

    //Distribución de energía
    for (int t = 0; t< DUR_SIM; t++) {// Instancia de tiempo
        for (int k = 0; k < numTri; k++) {// Triangulo 1
            if (MET.d[k][t]!=0){
                int tri, tim;
                //Energía de la Sala
                for (tri = 0; tri < numTri; tri++) {
                    if(matRoomTime.m[k][tri]!=0){
                        tim= matRoomTime.m[k][tri] + t;
                    }
                    if (tim < DUR_SIM) {
                        MET.d[tri][tim] += (MET.d[k][t] * matRoomAngles.d[k][tri]*(1-alfa));
                    }
                }

                //Energía del receptor
                for (int re = 0; re < r.NR; re++) {
                    if(matRecTime.m[re][k]!=0){
                        tim= matRecTime.m[re][k] + t;
                    }
                    if (tim < DUR_SIM) {
                        r.rec[re].eR[tim] +=  (MET.d[k][t] * matRecAngles.d[re][k]);
                    }
                };
            }
        }
    }
    /*
    eneRes = eneRay;
    for(int j=0;j<s.NRAYS;j++){
        //cout<<endl<<"Las energias residual es: "<<endl;
        for (int i = 1; i < reflexiones[j].N; i++) {
            //cout<<"Energia res "<<i<<": "<<eneRes<<endl;
            eneRes = eneRes * (1 - alfa) * (1 - delta); // nueva energia que se va a transmitir en la nueva reflexion
        }
    }*/
    cout << "*******************************************************\n"<<endl;

    cout<<"********Guardando matriz receptor en el archivo********"<<endl;
    ofstream archivoCSV1;
    guardar_archivo("datosReceptor.csv", archivoCSV1);
    for (int i = 0; i < r.NR; ++i) {
        for (int j = 0; j < DUR_SIM; ++j) {
            archivoCSV1 << r.rec[i].eR[j] << " ";
        }
        archivoCSV1 << endl;
    }
}
