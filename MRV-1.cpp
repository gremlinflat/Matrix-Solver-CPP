#ifndef MRV1
#define MRV1
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

void File_Counter(string importFile, int *baris, int *kolom, char tipe_matriks);

class Matrix{
    public :
    int m, n;
    ifstream myfile;
    char tipe_matrix;
    float temp=0;
    float dataMatrix[100][100];
    float Hasil[100];
    string namaFile;
    int flag;
    float det=1;
    int tukarCount=1;
    int tukarbaris=0;
    float tHasil[100];
    //constructor
    Matrix(int _m, int _n, char type){
                if(type == 'a'){
                        m = _m;
                        n = _n+1;
                        tipe_matrix = 'a';
                }else if(type == 'n'){
                        m = _m;
                        n = _n;
                        tipe_matrix = 'n';
                }
        }
    void inputfromTerminal(){
                for(int i=0;i<m;i++){
                        for(int j=0;j<n;j++){
                            
                                cin>>Matrix::dataMatrix[i][j];
                        }
                }
        }
    void inputfromFile(string importFile){
                //read baris & kolom
                
                myfile.open(importFile);
                for(int i=0;i<m;i++){
                        for(int j=0;j<n;j++){
                                myfile >> Matrix::dataMatrix[i][j];
                        }
                }
                myfile.close();
        }
     
    void printMatrix(string A){
        int stop, start=0;
        if(A == "invers"){
            stop = 2*n;
            start = n;
        }else if(A == "denganIdentitas"){
            start =0;
            stop = 2*n;
        }else{
            start = 0;
            stop = n;
        }
        cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
        for(int i=0;i<m;i++){
            for(int j=start;j<stop;j++){
                cout<<Matrix::dataMatrix[i][j]<<"\t";
            }
            cout<<"|"<<endl;
        }
        cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<endl;

    }
        
    void SPL_cetakLinear(){
                cout << "Persamaan Linear dari Matriks Eselon diatas adalah"<<endl;
                for (int i=0; i<m; i++){
                        int koef = 1;
                        for(int j=0;j<n;j++){
                                if(dataMatrix[i][j] != 0){
                                        if(j < n-1){
                                                if( (j !=0) && (dataMatrix[i][j]<0) ){
                                                        cout<<"(";
                                                }
                                                cout<<dataMatrix[i][j]<<"X"<<koef;
                                                if( (j !=0) && (dataMatrix[i][j]<0) ){
                                                        cout<<")";
                                                }
                                        }
                                        if( (j+1 < n-1) && dataMatrix[i][j+1]!=0){
                                                cout<<" + ";
                                        }else if(j+1 == n){
                                                cout<<" = "<<dataMatrix[i][j];
                                        }else{
                                                
                                        }
                                }
                                koef++;
                        }
                        cout<<endl;
                }
                cout<<endl;
        }

    bool SPL_CekSolusi(){
                bool a = true;
                for(int i=0;i<m;i++){
                        if( (abs(dataMatrix[i][i]) == 0) && (abs(dataMatrix[i][n-1]))){
                                cout<<"Tidak Memiliki solusi / Memiliki Solusi Trivial";
                                return false;
                        }
                }
                return a;
        }
    bool SPL_CekBarisNol(int a, int b){
        for(int i=0;i<n;i++){
            if( abs(dataMatrix[a][i]) == 0){
                return false;
            }
        }
    }
    bool Det_CekSolusi(){
                bool a = true;
                for(int i=0;i<m;i++){
                        if(dataMatrix[i][i] == 0){
                                cout<<"Matriks ini tidak memiliki determinan";
                                return false;
                        }
                }
                return a;
        }
    void SPL_PrintHasil(){
                cout<<"Solusi Dari Persamaan Linear ini adalah"<<endl;
                for (int i = 0; i < m; i++){
                        cout << "X[" << i+1 << "] = " << Hasil[i]<<endl;
                }
                cout<<endl;
        }

    void SPL_BackwardSubs(){
        for(int i=m-1;i>=0;i--){
            Hasil[i] = dataMatrix[i][n-1];
            for(int j=i+1; j<n; j++){
                Hasil[i] -=(dataMatrix[i][j] * Hasil[j]); //aljabar simpel
            }
                    //reassign variable hasil ke satuan
            Hasil[i] = Hasil[i] / dataMatrix[i][i];
        }
    }

    void SPL_DiagonaltoIdentity(){
                for(int i=0; i<m;i++){
                        float c = (1/dataMatrix[i][i]);
                        if(c != 1){
                                Matrix::OBE_kali(i, c);
                                cout<<"R"<<i+1<<" = R"<<i+1<<" / "<<1/c<<endl;
                                printMatrix("biasa");
                        }
                }
        }

    void SPL_GaussTukar(){
        for(int a=0;a<m;a++){ //Kolom
            int maxId = a;
            float maxValue = dataMatrix[a][a];
            for(int b=a+1;b<n-1;b++){ //Baris
                if( (dataMatrix[b][a]>0 ? dataMatrix[b][a] : -1*dataMatrix[b][a]) > maxValue){ //finding the biggest
                    maxValue = dataMatrix[b][a];
                    maxId = b;
                }
            }
            if(maxId != a){ // if tuker
                OBE_tukar(maxId, a);
                cout<<"TUKAR R"<<a+1<<" <-> R"<<maxId+1<<endl;
                printMatrix("biasa");
            tukarbaris++;
            }
        }
    }
    void SPL_GaussOperate(){
        for(int i = 0;i<m;i++){
            int mark = i;
            for(int j = i+1;j<n;j++){
                float c = (dataMatrix[j][i]/dataMatrix[mark][mark]);
                if(c !=0){
                    OBE_Operate(j, mark, c*-1 );
                    cout<<"R"<<j+1<<" = R"<<j+1<<" - "<<c<<"R"<<mark+1<<endl;
                    printMatrix("biasa");
                }
            }
        }
        // 
    }       
    void SPL_Jordan(){
        for(int i = m-1;i>0;i--){
            int mark = i;
            for(int j = i-1; j>=0; j--){
                if(i!=j){
                    float c = (dataMatrix[j][i]/dataMatrix[mark][mark]);
                    if(c != 0){
                        OBE_Operate(j, mark, c*-1 );
                        cout<<"R"<<j+1<<" = R"<<j+1<<" - "<<c<<"R"<<mark+1<<endl;
                        SPL_0Overflow("Matriks");
                        SPL_0Overflow("Hasil");
                        printMatrix("biasa");

                    }
                }       
            }
        }
        SPL_DiagonaltoIdentity();
        //assign hasil
        for (int i = 0; i < m; i++){
                Hasil[i] = dataMatrix[i][n-1];
        }            
    }
        
    void OBE_tukar(int R1, int R2){
                for (int i = 0; i < n; i++){
                        temp = dataMatrix[R1][i];
                        dataMatrix[R1][i] = dataMatrix[R2][i];
                        dataMatrix[R2][i] = temp;
                }
                
        }
    void OBE_Operate(int R1, int R2, float c){
                for (int i = 0; i < n; i++){
                        dataMatrix[R1][i] += (c*dataMatrix[R2][i]);                        
                }
                
        }
    void OBE_kali(int R1, float c){
                for (int i = 0; i < n; i++){
                        dataMatrix[R1][i] *=c;
                }
                
        }

    void DeterminanDiagonal(){//CARI DETERMINAN MATRIKS SEGITIGA ATAS/BAWAH
        float det=1;
        cout<<"Determinan Matriks ini adalah "<< endl <<"Det (A) = ";
        for(int i=0;i<m;i++){
            cout << dataMatrix[i][i]; 
            if(i<m-1){
                cout<<" * ";
            }
            det = det * dataMatrix[i][i];
        }
        det = pow(-1, tukarCount)*det;
        cout<<" * (-1)^"<<tukarCount;
        cout<<" = "<< det;
    }

    void Det_Cepat(){//mencari determinan dengan metode gauss untuk mendapatkan matriks segitiga tanpa mencetak ke layar <<BUAT FUNGSI CRAMMER>>
        det = 1;
        for(int a=0;a<m;a++){ //Kolom
            int maxId = a;
            float maxValue = dataMatrix[a][a];
            for(int b=a+1;b<n-1;b++){ //Baris
                if( (dataMatrix[b][a]>0 ? dataMatrix[b][a] : -1*dataMatrix[b][a]) > maxValue){
                    maxValue = dataMatrix[b][a];
                    maxId = b;
                }
            }
            if(maxId != a){
                OBE_tukar(maxId, a);
            }
        }

        for(int i = 0;i<m;i++){
            int mark = i;
            for(int j = i+1;j<n;j++){
                float c = (dataMatrix[j][i]/dataMatrix[mark][mark]);
                if(c !=0){
                    OBE_Operate(j, mark, c*-1 );
                }
            }
        }
        for(int i=0;i<m;i++){
            det = det * dataMatrix[i][i];
        }
        det = pow(-1, tukarCount) * det;
    }

    void inverse(){
        for(int i=0;i<m;i++){
            for(int j=n;j<2*n;j++){
                if(j-(n) == i){
                    dataMatrix[i][j]=1;
                }else{
                    dataMatrix[i][j]=0;
                }
            }
        }
        printMatrix("denganIdentitas");
        n = (2*n);
        SPL_GaussTukar();
        SPL_GaussOperate();
        SPL_Jordan();
        printMatrix("Invers");
        n = n/2;

    }

    void SPL_invers(){
        for (int i = 0; i < m; i++){ // perkalian matirks
            Hasil[i]=0;
            for (int j = 0; j < n; j++){
                Hasil[i] += dataMatrix[i][j] * dataMatrix[j][n-1];
            }
        }
    }

    void toNormal(){
        if(tipe_matrix != 'n'){
            n = m;
        }
    }
    void toAug(){
        if(tipe_matrix != 'a'){
            n = m;
        }
    }
    void grabtempHasil(){
        for(int i=0;i<m;i++){
            tHasil[i]=dataMatrix[i][n-1];
        }
    }
    void copyDataFrom(Matrix *data, int m2, int n2){//mengambil nilai matriks ke tempMatriks
        for(int i=0;i<m2;i++){
            for(int j=0;j<n2;j++){
                dataMatrix[i][j] = data->dataMatrix[i][j];
            }
        }
    }

    void SPL_0Overflow(string b){ //bug c++ kalo perkalian datamatriksnya lebih dari exponen -8
        if(b == "Matriks"){
            for(int i=0;i<m;i++){
                for(int j=0;j<n;j++){
                    if(abs(dataMatrix[i][j]) <= 0.1e-05){
                    //dataMatrix[i][j] = 0;
                    }
                }
            }
        }else if(b == "Hasil"){
            for(int j=0;j<n;j++){
                if(abs(Hasil[j]) <= 0.1e-05){
                    //Hasil[j] = 0;
                }
            }
        }else if(b == "tHasil"){
            for(int j=0;j<n;j++){
                if(abs(tHasil[j]) <= 0.1e-05){
                    //tHasil[j] = 0;
                }
            }
        }else if(b=="gauss"){
            for(int i=1;i<m;i++){
                for(int j=0;j<i;j++){
                    if(abs(dataMatrix[i][j]) <= 0.1e-05){
                        dataMatrix[i][j] = 0;
                    }
                }
            }
        }
    }
    void expansi_kofaktor(int ){
        det=1;
    }
    void SPL_Gauss_LastLeadingOne(){
        OBE_kali( (m-1), (1/dataMatrix[m-1][n-2]) );
    }

    

};
void inputMethod(char a, int *x, int *y, string *fName, bool *fromFile){
    cout<<endl<<"1. Input Matriks dari Keyboard "<<endl;
    cout      <<"2. Input Matriks dari file"<<endl;
    inputMet:
    cout<<endl<<"Pilihan <1-2> : ";
    int opsi11;
    cin>>opsi11;
    if(a == 'a'){
        switch(opsi11){
            case 1 :{//terminal
                cout<<"Masukan Jumlah Persamaan [m] = "; cin>>*x;
                cout<<"Masukan Jumlah Variable [n] = "; cin>>*y;
                *fromFile = false;
                break;
            }
            case 2 :{//file
                //string *fName = "matrix.txt";
                cout<<"Masukan Nama file (matrix.txt)= ";
                cin >> *fName;
                *fromFile = true;
                string namafiles = *fName;
                File_Counter(namafiles, x, y, a);
                break;
            }
            default :
                cout<<"Input tidak Valid"<<endl;
                goto inputMet;
        }
    }
    else if(a == 'n'){
        switch(opsi11){
            case 1 :{//terminal
                cout<<"Masukan ordo Matrix [n] = "; cin>>*x;
                *y = *x;
                *fromFile = false;
                break;
            }
            case 2 :{
                //string *fName = "matrix.txt";
                cout<<"Masukan Nama file (matrix.txt)= ";
                cin >> *fName;
                *fromFile = true;
                string namafiles = *fName;
                File_Counter(namafiles, x, y, a);
                break;
            }
            default :
                cout<<"Input tidak Valid"<<endl;
                goto inputMet;
        }
    }
}
void File_Counter(string importFile, int *baris, int *kolom, char tipe_matrix){
    ifstream file(importFile);
    string line;
    getline(file,line);         //baca string
    stringstream s;
    s << line;                   //masukin ke operasi string
    *kolom = 0;
    *baris = 1;
    float value;
    while (getline(file, line)){
            *baris = *baris+1; //setiap ketemu akhir string, tambah
    }
    while(s >> value){
            *kolom = *kolom+1;//ada sesuatu di baris selain spasi, tambahkan
    }
    if(tipe_matrix == 'a'){
        *kolom = *kolom-1;
    }
}

void Determinan_Menu(){
    cout<<"1. Metode Reduksi Baris"<<endl;
    cout<<"2. Metode Expansi Kofaktor"<<endl;
    int opsidet;
    opsidet:
    cout<<"Masukan Pilihan [1-2] = "; 
    cin>>opsidet;
    switch (opsidet){
        case 1:{ //obe reduksi
            int m1,n1;
            string FileName="matrix.txt";
            bool fromFile;
            char T='n';
            inputMethod(T, &m1, &n1, &FileName, &fromFile);
            Matrix Matriks(m1, n1, T);
            if(fromFile){
                Matriks.inputfromFile(FileName);
            }else{
                Matriks.inputfromTerminal();
            }
            if(m1!=n1){
                cout<<"Matriks harus berbentuk persegi"<<endl;
                return ;
            }
            Matriks.printMatrix("biasa");
            Matriks.SPL_GaussTukar();
            Matriks.SPL_GaussOperate();
            bool validitas = Matriks.Det_CekSolusi();
            if(validitas){
                Matriks.DeterminanDiagonal();
            }
        break;
        }
        case 2:{ //ekspansi kofaktor
            // //startinput
            // int m1,n1;
            // string FileName="matrix.txt";
            // bool fromFile;
            // char T='a'; //augmented
            // inputMethod(T, &m1, &n1, &FileName, &fromFile);
            // Matrix Matriks(m1, n1, T);
            // if(fromFile){
            //     Matriks.inputfromFile(FileName);
            // }else{
            //     Matriks.inputfromTerminal();
            // }
            // //endinput
            cout<<"MAAF PAK, UNTUK FUNGSI INI LOGIKA KAMI BELUM SAMPAI DAN BUTUH WAKTU LEBIH LAMA LAGI UNTUK MENYELESAIIN FUNGSI INI :(";
            break;
        }
        default :
            cout<<"input tidak valid";
            goto opsidet;
    }
}

void Invers_Menu(){
    //startinput
    int m1,n1;
    string FileName="matrix.txt";
    bool fromFile;
    char T='n'; //normal
    inputMethod(T, &m1, &n1, &FileName, &fromFile);
    Matrix Matriks(m1, n1, T);
    if(fromFile){
        Matriks.inputfromFile(FileName);
    }else{
        Matriks.inputfromTerminal();
    }
    //stopinput
    if(m1!=n1){
        cout<<"Matriks harus berbentuk persegi"<<endl;
        return ;
    }
    Matriks.inverse();
    cout<<endl<<"Matriks Balikan dari Matriks Tersebut adalah "<<endl;
    Matriks.printMatrix("invers");

}

void SPL_MENU(){
    cout<<endl<<"1. Metode Eliminasi Gauss"<<endl;
    cout<<"2. Metode Eliminasi Gauss-Jordan"<<endl;
    cout<<"3. Metode Matriks Balikan"<<endl;
    cout<<"4. Kaidah Crammer"<<endl;
    splmenu:
    int opsi;
    cout<<endl<<"Masukan Pilihan [1-4] = "; cin>>opsi;
    switch (opsi){  
        case 1 :{//gauss
            string FileName;
            int m1=0,n1=0;
            bool fromFile;
            char T='a'; //augmented
            inputMethod(T, &m1, &n1, &FileName,&fromFile); //dapetin baris kolom file
            Matrix Matriks(m1, n1, T);

            if(fromFile){
                Matriks.inputfromFile(FileName);
                Matriks.printMatrix("biasa");
            }else{
                Matriks.inputfromTerminal();
            }
            

            Matriks.SPL_GaussTukar();
            Matriks.SPL_GaussOperate();
            Matriks.SPL_0Overflow("gauss");
            bool validitas = Matriks.SPL_CekSolusi();
            if(validitas){
                    Matriks.SPL_cetakLinear();
                    Matriks.SPL_Gauss_LastLeadingOne();
                    Matriks.SPL_BackwardSubs();
                    Matriks.SPL_PrintHasil();
            }
            break;
        }   
        case 2 :{//gauss-jordan
            
            //startinput
            string FileName;
            int m1=0,n1=0;
            bool fromFile;
            char T='a'; //augmented
            inputMethod(T, &m1, &n1, &FileName,&fromFile); //dapetin baris kolom file
            Matrix Matriks(m1, n1, T);
            if(fromFile){
                Matriks.inputfromFile(FileName);
                Matriks.printMatrix("biasa");
            }else{
                Matriks.inputfromTerminal();
            }
            //stopinput

            Matriks.SPL_GaussTukar();
            Matriks.SPL_GaussOperate();
            
            bool validitas = Matriks.SPL_CekSolusi();
            if(validitas){
                Matriks.SPL_Gauss_LastLeadingOne();
                Matriks.SPL_Jordan();
                Matriks.SPL_0Overflow("Hasil");
                Matriks.SPL_PrintHasil();
            }
            break;
        }
        case 3 :{//Invers 
            //startinput
            string FileName;
            int m1=0,n1=0;
            bool fromFile;
            char T='a'; //augmented
            inputMethod(T, &m1, &n1, &FileName,&fromFile); //dapetin baris kolom file
            Matrix Matriks(m1, n1, T);
            if(fromFile){
                Matriks.inputfromFile(FileName);
                Matriks.printMatrix("biasa");
            }else{
                Matriks.inputfromTerminal();
            }
            //stopinput

            Matriks.grabtempHasil();
            Matriks.toNormal();

            float tempInvers[m1][n1];

            Matriks.inverse();
            int stop = ((n1)*2);
            int start = n1;
            int counter1=0;
            int counter2=0;
            
            for(int i=0;i<m1;i++){ //mengambil matriks balikan dan memasukan ke matriks temporary
                counter2=0;
                for(int j=start;j<n1*2;j++){
                    tempInvers[counter1][counter2] = Matriks.dataMatrix[i][j];
                    counter2++;
                }
                counter1++;
            }

            Matriks.printMatrix("invers");
            cout<<endl<<"Maka, "<<endl;

            for(int i=0;i<m1;i++){//loop ngeprint
                for(int j=0;j<n1;j++){
                    cout<<tempInvers[i][j]<<"\t";
                }
                if( i==((m1)/2) ){
                    cout<<"   "<< "X  ";
                }
                cout<<"\t \t"<<Matriks.tHasil[i]<<endl;
            }

            for (int i = 0; i < m1; i++){//perkalian matriks
                Matriks.Hasil[i] = 0;
                for (int j = 0; j < n1; j++){
                    Matriks.Hasil[i] += tempInvers[i][j] * Matriks.tHasil[j];
                }       
            }
            cout<<endl;
            Matriks.SPL_PrintHasil();
            break;
            }
        case 4 : {//crammer
            //startinput
            string FileName;
            int m1=0,n1=0;
            bool fromFile;
            char T='a'; //augmented
            inputMethod(T, &m1, &n1, &FileName,&fromFile); //dapetin baris kolom file
            Matrix Matriks(m1, n1, T);
            if(fromFile){
                Matriks.inputfromFile(FileName);
                Matriks.printMatrix("biasa");
            }else{
                Matriks.inputfromTerminal();
            }
            //stopinput

            Matriks.grabtempHasil();
            Matriks.toNormal();

            Matrix D1(m1, n1, T);
            D1.toNormal();

            D1.copyDataFrom(&Matriks, m1, n1);
            
            float subDeterminan[n1];
            
            D1.Det_Cepat();
            float DeterminanUtama = D1.det;
            
            if(abs(D1.det) == 0){
                cout<<"TIDAK MEMILIKI Solusi Karena Determinan Utamanya Adalah Nol";
                return;
            }
            
            Matrix tempMatriks(m1, n1, T);
            tempMatriks.toNormal();

            for(int i=0;i<m1;i++){
                tempMatriks.copyDataFrom(&Matriks, m1, n1);//reset data temp;
                for(int j=0;j<n1;j++){//tukar kolom
                    tempMatriks.dataMatrix[j][i] = Matriks.tHasil[j];
                }
                cout<<"A"<<i+1<<endl;
                tempMatriks.printMatrix("biasa");
                tempMatriks.Det_Cepat();
                subDeterminan[i]=tempMatriks.det;
                //cout<<"Det |A"<<i+1<<"| = "<<subDeterminan[i]<<endl;
            }
            cout<<"Det |A| = "<<DeterminanUtama<<endl;
            for(int i=0;i<m1;i++){
                cout<<"Det |A"<<i+1<<"| = "<<subDeterminan[i]<<" \t";
                if(i%3==2){
                    cout<<endl;
                }
                Matriks.Hasil[i]=subDeterminan[i]/DeterminanUtama;
            }
            cout<<"Solusi Dari Persamaan Linear ini adalah"<<endl;
            for (int i = 0; i < m1; i++){
                cout << "X[" << i+1 << "] = " <<subDeterminan[i]<<"/"<<DeterminanUtama<<" = " <<Matriks.Hasil[i]<<endl;
            }
            cout<<endl;
            break;
            }
        default :
                cout<<"Input tidak Valid"<<endl;
                goto splmenu;
    }
}



#endif