/*****************************************************************************************************************************
 * Proyecto    : Física Computacional II.
 * Version     : 17/07/2023 (running ok).
 * Descripcion : Ley de enfriamiento de Newton (simulación t = -ln[(T-Ta)/(To-Ta)]/k).
 * Autor       : Aros D., Campaña B., Jurado Ordoñez Y., Palacios A.
 *****************************************************************************************************************************/

using namespace std;
using namespace TMath;

// Parámetros globales.
Double_t Ta = -8.;                            // Temperatura ambiente [ºC].
Double_t To = 74;                             // Temperatura agua [ºC].
Double_t k  = 0.000764;                       // Constante de enfriamiento.

// Constructores (sirven para inicializar un objeto y establecer sus propiedades y valores predeterminados).
Double_t len_dif(Double_t x, Double_t y);
Double_t rk4_solver(Double_t xo, Double_t yo, Double_t h, Double_t x, Double_t len_dif(Double_t xx, Double_t yy));
Double_t fitFunc(Double_t* x, Double_t* par);
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny, Float_t lMargin, Float_t rMargin,Float_t bMargin, Float_t tMargin);

// Función principal:
void gr2(){
    // Información del experimento ...........................................................................................
    const Int_t npts  	= 10;                         // Número de puntos para las graficas.
  
    // Declaración de variables ..............................................................................................
    
    Double_t tiempo_plas_hab[10]      = {0., 61.012, 154.066, 271.057, 426.037, 635.080, 880.026, 1227.030, 1693.073, 2451.036};
    Double_t tiempo_plas_nev[10]      = {0., 49.046, 123.021, 214.028, 322.074, 591.064, 762.054, 981.024, 1241.062, 1541.059};
    Double_t tiempo_porc_hab[10]      = {0., 31.083, 105.037, 193.050, 325.006, 498.057, 734.031, 1040.016, 1478.095, 2105.031};
    Double_t tiempo_porc_nev[10]      = {0., 60, 155.051, 268.006, 334.0.38, 470.024, 577.038, 768.061, 953.062, 1185.045};
    Double_t tiempo_vidr_hab[10]      = {0., 55.041, 140.026, 253.067, 399.059, 582.015, 830.067, 1157.062, 1622.082, 2318.051};
    Double_t tiempo_vidr_nev[10]      = {0., 48.069, 116.010, 201.013, 299.082, 415.034, 545.016, 699.061, 896.004, 1145.097};
    Double_t temperatura_real[10]     = {74., 70., 65., 60., 55., 50., 45., 40., 35., 30.};
    Double_t temperatura_real_p[10]   = {71, 66, 59, 53., 50., 45., 40., 35., 30.,25};
    Double_t tiempo_real_err[10]      = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
    Double_t temperatura_real_err[10] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
    
    // Cálculos ..............................................................................................................
    
    TF1 *f1 = new TF1("f1",fitFunc, 0., 3000.,2);
    f1->SetParNames("To","k");
    f1->SetParameters(To,k);
    
    // Graficas .........................................................................................................
    
    TCanvas *C = (TCanvas*) gROOT->FindObject("C");
    if (C) delete C;
    C = new TCanvas("C","canvas",1024,640);
    C->SetFillStyle(4000);
    
    // Numero de PADS
    const Int_t Nx = 2;
    const Int_t Ny = 2;
    
    // Mergenes
    Float_t lMargin = 0.1409;
    Float_t rMargin = 0.1209;
    Float_t bMargin = 0.08;
    Float_t tMargin = 0.05;
    
    // Canvas setup
    CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
    
    TPad *pad[Nx][Ny];
    
    for (Int_t i=0;i<Nx;i++) {
        for (Int_t j=0;j<Ny;j++) {
            C->cd(0);
            
            // Get the pads previously created.
            char pname[16];
            sprintf(pname,"pad_%i_%i",i,j);
            pad[i][j] = (TPad*) gROOT->FindObject(pname);
            pad[i][j]->Draw();
            pad[i][j]->SetFillStyle(4000);
            pad[i][j]->SetFrameFillStyle(4000);
            
        }
    }
    
    Float_t xFactor = 0.;
    Float_t yFactor = 0.;
    
    ///////////////////////////////////////////////////////////
    pad[0][1]->cd();                     // Crea el 1er pad.
    
    xFactor = pad[0][0]->GetAbsWNDC()/pad[0][1]->GetAbsWNDC();
    yFactor = pad[0][0]->GetAbsHNDC()/pad[0][1]->GetAbsHNDC();
    
    TGraphErrors *gr1 = new TGraphErrors(npts, tiempo_plas_hab, temperatura_real, tiempo_real_err, temperatura_real_err);
    gr1->SetMarkerColor(kBlue);
    gr1->SetMarkerStyle(8);
    gr1->SetMarkerSize(1);
    gr1->Draw("ap");
    
    TGraphErrors *gr2 = new TGraphErrors(10, tiempo_plas_nev, temperatura_real, tiempo_real_err, temperatura_real_err);
    gr2->SetMarkerColor(kGreen+3);
    gr2->SetMarkerStyle(8);
    gr2->SetMarkerSize(1);
	gr2->Fit("f1");
    gr2->Draw("p");    
    
    gr1->GetXaxis()->SetRangeUser(0,3000);
    //gr1->GetYaxis()->SetRangeUser(0,2000);
    gr1->GetXaxis()->SetTitle("Tiempo [s]");
    gr1->GetYaxis()->SetTitle("Temperatura [^{o}C]");
	gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
    
    TLegend* leg1 = new TLegend(0.5, 0.6, 0.85, 0.8);
    leg1->SetBorderSize(0);
    leg1->SetHeader("   Agua-Plastico");
    leg1->AddEntry(gr1,"Datos Habitacion ","pe");
    leg1->AddEntry(gr2,"Datos Nevera","pe");
    leg1->Draw();
    
    
    
    /////////////////////////////////////////////////////////
    pad[1][1]->cd();                   // Crea el 2do pad.
    
    xFactor = pad[0][0]->GetAbsWNDC()/pad[1][1]->GetAbsWNDC();
    yFactor = pad[0][0]->GetAbsHNDC()/pad[1][1]->GetAbsHNDC();
    
    TGraphErrors *gr3 = new TGraphErrors(npts, tiempo_porc_hab, temperatura_real, tiempo_real_err, temperatura_real_err);
    gr3->SetMarkerColor(kBlue);
    gr3->SetMarkerStyle(8);
    gr3->SetMarkerSize(1);
    gr3->Draw("ap");
    
    TGraphErrors *gr4 = new TGraphErrors(10, tiempo_porc_nev, temperatura_real_p, tiempo_real_err, temperatura_real_err);
    gr4->SetMarkerColor(kGreen+3);
    gr4->SetMarkerStyle(8);
    gr4->SetMarkerSize(1);
    gr4->Draw("p");
	gr4->Fit("f1");
	
    gr3->GetXaxis()->SetRangeUser(0,3000);
    //g3->GetYaxis()->SetRangeUser(0,2000);
    gr3->GetXaxis()->SetTitle("Tiempo [s]");
    gr3->GetYaxis()->SetTitle("Temperatura [^{o}C]");
	gr3->GetXaxis()->CenterTitle();
    gr3->GetYaxis()->CenterTitle();
    
    TLegend* leg2 = new TLegend(0.5, 0.6, 0.85, 0.8);
    leg2->SetBorderSize(0);
    leg2->SetHeader("Agua-Porcelana");
    leg2->AddEntry(gr3,"Datos Habitacion","pe");
    leg2->AddEntry(gr4,"Datos Nevera","pe");
    leg2->Draw();
	
    
    ///////////////////////////////////////////////////////////
    pad[0][0]->cd();                     // Crea el 3er pad.
    
    xFactor = pad[0][0]->GetAbsWNDC()/pad[0][0]->GetAbsWNDC();
    yFactor = pad[0][0]->GetAbsHNDC()/pad[0][0]->GetAbsHNDC();
    
	TGraphErrors *gr5 = new TGraphErrors(npts, tiempo_vidr_hab, temperatura_real, tiempo_real_err, temperatura_real_err);
    gr5->SetMarkerColor(kBlue);
    gr5->SetMarkerStyle(8);
    gr5->SetMarkerSize(1);
    gr5->Draw("ap");
    
    TGraphErrors *gr6 = new TGraphErrors(10, tiempo_vidr_nev, temperatura_real, tiempo_real_err, temperatura_real_err);
    gr6->SetMarkerColor(kGreen+3);
    gr6->SetMarkerStyle(8);
    gr6->SetMarkerSize(1);
    gr6->Draw("p");
	gr6->Fit("f1");
	
    gr5->GetXaxis()->SetRangeUser(0,3000);
    //gr5->GetYaxis()->SetRangeUser(0,2000);
    gr5->GetXaxis()->SetTitle("Tiempo [s]");
    gr5->GetYaxis()->SetTitle("Temperatura [^{o}C]");
	gr5->GetXaxis()->CenterTitle();
    gr5->GetYaxis()->CenterTitle();
    
    TLegend* leg2 = new TLegend(0.5, 0.6, 0.85, 0.8);
    leg2->SetBorderSize(0);
    leg2->SetHeader("Agua-Vidrio");
    leg2->AddEntry(gr5,"Datos Habitacion","pe");
    leg2->AddEntry(gr6,"Datos Nevera","pe");
    leg2->Draw();

    	
    ////////////////////////////////////////////////////////////
    pad[1][0]->cd();                      // Crea el 4to pad.
    
    xFactor = pad[0][0]->GetAbsWNDC()/pad[1][0]->GetAbsWNDC();
	yFactor = pad[0][0]->GetAbsHNDC()/pad[1][0]->GetAbsHNDC();
		
	TF1 *f2 = gr2->GetFunction("f1");
	f2->SetLineColor(kCyan-3);
	f2->Draw();
		
	TF1 *f3 = gr4->GetFunction("f1"); // Corrección aquí para obtener la len_difción de ajuste de gr4
	f3->SetLineColor(kOrange-2);
	f3->Draw("same");
		
	TF1 *f4 = gr6->GetFunction("f1"); // Corrección aquí para obtener la len_difción de ajuste de gr6
	f4->SetLineColor(kGray+2);    
	f4->Draw("same");
		
	f2->GetXaxis()->SetRangeUser(-20,1250);
	// f2->GetYaxis()->SetRangeUser(0,2000);
	f2->GetXaxis()->SetTickLength(0.06);
	f2->GetXaxis()->SetTitle("Tiempo [s]");
	f2->GetYaxis()->SetTitle("Temperatura [ ^{o}C ]");
	f2->GetXaxis()->CenterTitle();
	f2->GetYaxis()->CenterTitle();
			
	TLegend* leg4 = new TLegend(0.5, 0.6, 0.85, 0.8);
	leg4->SetBorderSize(0);
	leg4->SetHeader("Mejores Ajustes");
	leg4->AddEntry(f2,"Agua-Plastico Nevera","l");
	leg4->AddEntry(f3,"Agua-Porcelana Nevera","l");
	leg4->AddEntry(f4,"Agua-Vidrio Nevera","l");
	leg4->Draw();
		
		
	C->cd(0);
	C->Update();
	C->Modified();
}

///////////////////////////////////////////   Funciones para el ajuste   ///////////////////////////////////////////

Double_t len_dif(Double_t x, Double_t y){
    return -k * (y - Ta);
}


Double_t rk4_solver(Double_t xo, Double_t yo, Double_t h, Double_t x){
    Int_t n = (x - xo)/h;
    Double_t y = yo;
		for (Int_t i = 0; i<n; i++) {
            Double_t k1 = h*len_dif(xo, y);
            Double_t k2 = h*len_dif(xo + 0.5*h, y + 0.5*k1);
            Double_t k3 = h*len_dif(xo + 0.5*h, y + 0.5*k2);
            Double_t k4 = h*len_dif(xo + h, y + k3);
            y += (k1 + 2.*k2 + 2.*k3 + k4)/6.;
            xo += h;
        }
    return y;
}	

Double_t fitFunc(Double_t* x, Double_t* par){
    Double_t To = par[0];
    Double_t k = par[1];

    Double_t t = x[0];
    Double_t y = rk4_solver(0, To, 0.1, t); // Calculamos el valor de y con Runge-Kutta

    return y;
}


////////////////////////////////////////////    Divición del canvas    //////////////////////////////////////////
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny, Float_t lMargin, Float_t rMargin,Float_t bMargin, Float_t tMargin){
    if (!C) return;
    
    // Setup Pad layout:
    Float_t vSpacing = 0.03;
    Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
    
    Float_t hSpacing = 0.03;
    Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
    
    Float_t vposd,vposu,vmard,vmaru,vfactor;
    Float_t hposl,hposr,hmarl,hmarr,hfactor;
    
    for (Int_t i=0;i<Nx;i++) {
        
        if (i==0) {
            hposl = 0.0;
            hposr = lMargin + hStep;
            hfactor = hposr-hposl;
            hmarl = lMargin / hfactor;
            hmarr = 0.0;
        } else if (i == Nx-1) {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep + rMargin;
            hfactor = hposr-hposl;
            hmarl = 0.0;
            hmarr = rMargin/hfactor;
        } else {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
            hfactor = hposr-hposl;
            hmarl = 0.0;
            hmarr = 0.0;
        }
        
        for (Int_t j=0;j<Ny;j++) {
            
            if (j==0) {
                vposd = 0.0;
                vposu = bMargin + vStep;
                vfactor = vposu-vposd;
                vmard = bMargin / vfactor;
                vmaru = 0.0;
            } else if (j == Ny-1) {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep + tMargin;
                vfactor = vposu-vposd;
                vmard = 0.0;
                vmaru = tMargin / vfactor;
            } else {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep;
                vfactor = vposu-vposd;
                vmard = 0.0;
                vmaru = 0.0;
            }
            
            C->cd(0);
            
            char name[16];
            sprintf(name,"pad_%i_%i",i,j);
            TPad *pad = (TPad*) gROOT->FindObject(name);
            if (pad) delete pad;
            pad = new TPad(name,"",hposl,vposd,hposr,vposu);
            pad->SetLeftMargin(hmarl);
            pad->SetRightMargin(hmarr);
            pad->SetBottomMargin(vmard);
            pad->SetTopMargin(vmaru);
            
            pad->SetFrameBorderMode(0);
            pad->SetBorderMode(0);
            pad->SetBorderSize(10);
            
            pad->Draw();
        }
    }
}
