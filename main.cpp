#include <iostream>
#include "MRV-1.cpp"
using namespace std;
/*
    Tugas Besar MRV 1

    Fahri Novaldi               119140205
    Ilham Nofri Yandra          119140133
    M. Ammar Fadhila Ramadhan   119140029



*/
int main(){
    cout<<"============ Menu ==========="<<endl;
    cout<<"= 1. Sistem Persamaan Linier="<<endl;
    cout<<"= 2. Determinan             ="<<endl;
    cout<<"= 3. Matriks Balikan        ="<<endl;
    cout<<"= 4. Keluar                 ="<<endl;
    cout<<"-----------------------------"<<endl;
    pilihan:
    cout<<"Pilihan <1-4> : ";int pilihan;cin>>pilihan;
    switch (pilihan) {
        
            case 1:
                SPL_MENU();
                break;
            case 2:
                Determinan_Menu();
                break;
            case 3:
                Invers_Menu();
                break;
            case 4:
                cout<<"Program Ditutup";
                break;
            default:
                cout<< "Opsi tidak valid"<<endl;
                goto pilihan;
                break;
        }
    return 0;
    
}