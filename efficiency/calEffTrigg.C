#include <TPRegexp.h>
#include <TSystem.h>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <AliJRunTable.cxx>
#include <vector>
#include <iostream>
#include <TNamed.h>

using namespace std;

enum { kJFilterTPCOnly, kJFilterRaa, kJFilterTightDCA, kJFilterGlobalDCA, kJFilterGlobalSDD, kJFilterHybrid,
    kJNFilter
};
TString strTrk[] = { "TPCOnly","Raa","Global+TightDCA","Global+DCA", "Global+SDD", "Hybrid",
};
enum { kJPhysicsPrimary, kJFake, kNPrimaryStatus };
enum { kJMCTrack, kJGlobal, kJTPCOnly, kJGCG , kNTrackType };
enum { kNVtxBin=2 };
enum { kNCentBin=22};
enum { kRE, kMC };
enum { kPP, kAA, kPA };

TH1D* EffRatio(TH1D *h0,TH1D* H1, TH1D* H2){
    int nb =H1->GetNbinsX();
    if( h0 != H1 && h0!=H2){
        h0->Reset();
    }

    for(int i=1;i<=nb;i++){
        double c1 = H1->GetBinContent(i);
        double c2 = H2->GetBinContent(i);
        double e1 = H1->GetBinError(i);

        if(c2!=0){
            h0->SetBinContent(i,c1/c2);
            h0->SetBinError(i,e1/c2);
        }else{
            h0->SetBinContent(i,0);
            h0->SetBinError(i,0);
        }
    }

    return h0;
}               

class JCorranEff {
    public:
        JCorranEff( TString filename ); //TODO
        void LoadHist();
        void Evaluate( int ivtx, int icent0, int icent1, int ifilter );
        void EvaluateAll();
        void SetCentBinOut(TString centbinstr);
        void Write();
        void WriteName(TString name, TString title, int i);
        TH1D * LoadRebinScale( const char * name );

        TFile * fInRoot;
        TString fRootName;
        TDirectory * fEffD;
        TFile * fOutRoot;
        TString fOutRootName;
        TDirectory * fEffDOut;
        int     fNVtxBin;
        int     fNCentBin;
        vector<double>  fCentBinOut;
        int     fNCentBinOut;
        int     fRebin;
        double *     fRebins;
        int     fBeamType;

        AliJRunTable * fRT;


        TH1D * hChargedPtMC[kNVtxBin][kNCentBin];
        TH1D * hChargedPtMCTrigg[kNVtxBin][kNCentBin];
        TH1D * hChargedPtMCTriggVtx[kNVtxBin][kNCentBin];
        TH1D * hChargedPtMCRecoCentVtx[kNVtxBin][kNCentBin][kJNFilter][kNPrimaryStatus][kNTrackType];
        TH1D * hChargedPtDataCentVtx[kNVtxBin][kNCentBin][kJNFilter][kNTrackType];
        TH2D * h2VtxCent;

        TH1D * hCor;
        TH1D * hEff;
        TH1D * hCon;

        TH1D * hEffTrigg;
        TH1D * hEffTriggVtx;

        TGraphErrors * gCor[kNCentBin][kJNFilter];
        TGraphErrors * gEff[kNCentBin][kJNFilter];
        TGraphErrors * gCon[kNCentBin][kJNFilter];
        TGraphErrors * gEffTrigg[kNCentBin][kJNFilter];
        TGraphErrors * gEffTriggVtx[kNCentBin][kJNFilter];
};

JCorranEff::JCorranEff( TString filename ):fRootName(filename){
    fNVtxBin = 1;
    fNCentBin = 1;
}

TH1D * JCorranEff::LoadRebinScale( const char * name ){
    TH1D *h = (TH1D*) fEffD->Get( name );
    if(!h){
        cout<<"J_ERROR : No "<<name<<endl;
        gSystem->Exit(1);
    }
    h->Print();
    h=(TH1D*)h->Rebin(fRebin,TString(h->GetName())+"_rebin", fRebins);
    h->Scale(1,"width");
    return h;
}

void JCorranEff::LoadHist(){
    cout<<"J_LOG :: Load Eff Hist"<<endl;
    if( fRT->IsPP() ) fNCentBin = 1;
    if( fRT->IsPA() ) fNCentBin = 1; //TEMP
    //fNCentBin = 1;
    fOutRoot = new TFile(fOutRootName,"RECREATE");
    fEffDOut = fOutRoot->mkdir("Efficiency");
    fInRoot = new TFile(fRootName,"READ");
    fEffD   = (TDirectory*) fInRoot->Get("JCORRANEfficiency");//EffHist");
    for( int ivtx=0;ivtx<fNVtxBin;ivtx++ ){
        for( int icent=0;icent<fNCentBin ;icent++ ){
            hChargedPtMCTriggVtx[ivtx][icent] = LoadRebinScale( Form("hChargedPtMCTriggVtx%02d%02d", ivtx, icent ) );
            hChargedPtMCTrigg[ivtx][icent] = LoadRebinScale( Form("hChargedPtMCTrigg%02d%02d", ivtx, icent ) );
            hChargedPtMC[ivtx][icent] = LoadRebinScale( Form("hChargedPtMC%02d%02d", ivtx, icent ) );
            for( int ifilter=0;ifilter<kJNFilter;ifilter++ ){
                for( int itt=0;itt<kNTrackType;itt++ ){
                    for( int ipri=0;ipri<kNPrimaryStatus;ipri++ ){
                       hChargedPtMCRecoCentVtx[ivtx][icent][ifilter][ipri][itt] = LoadRebinScale(  Form("hChargedPtMCRecoCentVtx%02d%02d%02d%02d%02d", ivtx, icent, ifilter, ipri, itt ) );
                    }
                }
            }
        }
    }
    cout << "Nvtx= " << fNVtxBin << ", NCentBin = " << fNCentBin << endl;
}

void PrintAxis( TAxis * ax ){
   cout<<"# of bins = "<<ax->GetNbins()<<endl; 
}

void JCorranEff::Evaluate( int ivtx, int icent0, int icent1, int ifilter ){
    cout<<"J_LOG : Evaluate "<<Form("ivtx:%d icent0:%d icent1:%d ifilter:%d ", ivtx, icent0, icent1, ifilter)<<strTrk[ifilter]<<endl;
    int itt=kJGlobal;
    if( ifilter == kJFilterTPCOnly ) itt=kJTPCOnly;

    TH1D * hMCTrue = (TH1D*)(hChargedPtMCTrigg[0][0]->Clone());
    hMCTrue->Reset();
    TH1D * hMCTrueTrigg = (TH1D*)(hChargedPtMCTrigg[0][0]->Clone());
    hMCTrueTrigg->Reset();
    TH1D * hMCTrueTriggVtx = (TH1D*)(hChargedPtMCTrigg[0][0]->Clone());
    hMCTrueTriggVtx->Reset();
    TH1D * hMCRecoPrim = (TH1D*)(hMCTrue->Clone());
    hMCRecoPrim->Reset();
    TH1D * hMCRecoFake = (TH1D*)(hMCTrue->Clone());
    hMCRecoFake->Reset();
    //for( int icent=icent0;icent<=icent1;icent++ ) // TODO ?
    for( int icent=icent0;icent<icent1;icent++ )
    {
        if ( icent >= kNCentBin ) continue; 
        cout << "Adding 2 hMCTrue, vtx=" << ivtx << " icent=" << icent << endl;
        hMCTrue->Add( hChargedPtMC[ivtx][icent] );
        cout << "Adding 2 hMCTrueTrig" << endl;
        hMCTrueTrigg->Add( hChargedPtMCTrigg[ivtx][icent] );
        cout << "Adding 2 hMCTrueTrigVtx" << endl;
        hMCTrueTriggVtx->Add( hChargedPtMCTriggVtx[ivtx][icent] );
        hMCRecoPrim->Add(hChargedPtMCRecoCentVtx[ivtx][icent][ifilter][kJPhysicsPrimary][kJMCTrack]);
        hMCRecoFake->Add(hChargedPtMCRecoCentVtx[ivtx][icent][ifilter][kJFake][itt]);

        //==== TODO Stat Error 
        //==== Efficiency = Np/Na
        hEff = (TH1D*) hMCRecoPrim->Clone(Form("hEff%02d%02d%02d",ivtx,icent,ifilter));
        EffRatio( hEff, hMCRecoPrim, hMCTrueTriggVtx ); 
        //==== Contamination = Nf/(Np+Nf) 
        hCon = (TH1D*) hMCRecoPrim->Clone(Form("hCon%02d%02d%02d",ivtx,icent,ifilter));
        hCon->Add( hMCRecoFake );
        EffRatio( hCon, hMCRecoFake, hCon );
        //==== Correction = (Np+Nf)/Na = Eff/(1-Con) = ( Np/Na * (Np+Nf)/Np)
        hCor = (TH1D*) hMCRecoPrim->Clone( Form("hCor%02d%02d%02d",ivtx, icent, ifilter));
        hCor->Add( hMCRecoFake );
        EffRatio( hCor, hCor, hMCTrueTriggVtx );

        //==== Efficiency = Np/Na
        hEffTrigg = (TH1D*) hMCTrueTrigg->Clone(Form("hEffTrigg%02d%02d%02d",ivtx,icent,ifilter));
        EffRatio( hEffTrigg, hEffTrigg, hMCTrue ); 
        hEffTriggVtx = (TH1D*) hMCTrueTriggVtx->Clone(Form("hEffTriggVtx%02d%02d%02d",ivtx,icent,ifilter));
        EffRatio( hEffTriggVtx, hEffTriggVtx, hMCTrueTrigg ); 

        cout << "Printing bin contents..." << endl;
        for( int i=0;i<hEff->GetNbinsX();i++ ){
            cout<<Form("%.3f", hEff->GetBinCenter(i+1)) << endl;//" " << Form("%-7.3f",hEff->GetBinContent(i-1));
        }
        cout<<endl;


        fEffDOut->cd();
        hEff->Write();
        hCor->Write();
        hCon->Write();
        hEffTrigg->Write();
        hEffTriggVtx->Write();
        gROOT->cd();
    }

}

void JCorranEff::EvaluateAll(){
    cout<<"J_LOG : Evaluate All"<<endl;
    int ivtx = 0;
    for( int iCentBinOut=0; iCentBinOut<fNCentBinOut;iCentBinOut++ )
    {
        double cent0 = fCentBinOut[iCentBinOut];
        double cent1 = fCentBinOut[iCentBinOut+1];

        //==== Get input bins from output bins
        int icent0 = 0;
        int icent1 = 0;
        icent0 = int(cent0/100.*fNCentBin);// 5% bins
        icent1 = int(cent1/100.*fNCentBin);
        for( int i=0;i<kJNFilter;i++ ){ //==== TRACK CUT LOOP
            fEffD->cd();
            cout<<"CENT BIN = "<<icent0<<"\t"<<icent1<<"\t"<<iCentBinOut<<endl;
            Evaluate(ivtx, icent0, icent1, i);
            hCor->Print();
            gCor[iCentBinOut][i] = new TGraphErrors(hCor);
            gEff[iCentBinOut][i] = new TGraphErrors(hEff);
            gCon[iCentBinOut][i] = new TGraphErrors(hCon);
            gEffTrigg[iCentBinOut][i] = new TGraphErrors(hEffTrigg);
            gEffTriggVtx[iCentBinOut][i] = new TGraphErrors(hEffTriggVtx);
        }
    }
}

void JCorranEff::SetCentBinOut(TString centbinstr){
    fCentBinOut.clear();
    TPMERegexp xreg("\\s+");
    xreg.Split( centbinstr );
    vector<double>centbinout;
    for( int i=0;i<xreg.NMatches();i++ ){
        TString nums = xreg[i]; 
        if( nums.Length() == 0 ) continue;
        double cent = nums.Atof();
        fCentBinOut.push_back( cent );
    }
    fNCentBinOut = fCentBinOut.size() - 1;
}

void JCorranEff::Write(){
    fEffDOut->cd();
    TAxis * xCentBin = new TAxis( fNCentBinOut, &(fCentBinOut[0]) );
    xCentBin->Write("CentralityBin");
    TAxis * xVtxBin = new TAxis( 1, -10, 10 );
    xVtxBin->Write("VtxBin");
    TAxis * xPtBin = hCon->GetXaxis();
    xPtBin->Write("PtBin");
    TAxis * xTrackCut = new TAxis( kJNFilter, 0, kJNFilter );
    WriteName( "BeamType", fRT->GetBeamTypeStr(), fRT->GetBeamType() );
    WriteName( "Period",fRT->GetPeriodName(), fRT->GetPeriod() );
    WriteName( "RunNumber",Form("%d",fRT->GetRunNumber()), fRT->GetRunNumber()) ;
    for( int i=0;i<kJNFilter;i++ ){
        xTrackCut->SetBinLabel(i+1, strTrk[i]);
    }
    xTrackCut->Write("TrackCut");

    for( int icent=0;icent<fNCentBinOut;icent++){
        for( int icut=0;icut<kJNFilter;icut++ ){
            gEff[icent][icut]->Write(Form("gEff%02d%02d%02d",0,icent,icut));
            gCon[icent][icut]->Write(Form("gCon%02d%02d%02d",0,icent,icut));
            gCor[icent][icut]->Write(Form("gCor%02d%02d%02d",0,icent,icut));
            gEffTrigg[icent][icut]->Write(Form("gEffTrigg%02d%02d%02d",0,icent,icut));
            gEffTriggVtx[icent][icut]->Write(Form("gEffTriggVtx%02d%02d%02d",0,icent,icut));

        }
    }
    fOutRoot->Write();
}

void JCorranEff::WriteName(TString name, TString title, int i){
    TNamed * named = new TNamed( name, title);
    named->SetUniqueID(i);
    named->Write();
}


void calEffTrigg(TString input,TString tag="" ){

    gROOT->LoadMacro("AliJRunTable.cxx+g");
    AliJRunTable rt;
    rt.ParseString(input);
    //TODO rt.Print();

    int runnumber       = rt.GetRunNumber();
    //rt.SetRunNumber( runnumber );
    //int period          = rt.GetPeriod();
    TString periodStr   = rt.GetPeriodName();
    TString periodMCStr = rt.GetPeriodMCName();
    int beamType        = rt.GetBeamType();
    /*
       AliJEfficiency effC;
       effC.SetDataPath("eff");
       effC.SetRunNumber(runnumber);
       effC.SetPeriod( periodStr );
       effC.SetMCPeriod( periodMCStr );
       effC.SetTag( tag );
       TString OutName = effC.GetEffFullName();
       */
    // Eff--LHC10b-LHC10d1-0-r1449.root
    TString OutName = Form("eff/Eff--%s-%s-%d-%s.root", periodStr.Data(), periodMCStr.Data(),runnumber, tag.Data());
    cout << OutName << endl;
    JCorranEff eff(input);
    eff.fOutRootName = OutName;
    eff.fBeamType = beamType;
    double xReb[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.5,1.8,2.1,2.4,2.7,3,3.3,4.5, 6.5,15};
    //double xReb[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.5,1.8,2.1,2.4,2.7,3,3.3,4.5, 6.5,30};
    //double xReb[]={0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3, 3.3, 4.5, 5.7, 7, 10, 15, 20, 30, 50};
    //    double xReb[]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.5,1.8,2.1,2.4,2.7,3,3.3,4.5,5.7,7,20};
    int nReb = sizeof(xReb)/sizeof(double)-1;
    eff.fRebin = nReb;
    eff.fRebins = xReb;
    if( rt.IsPP() || rt.IsPA() ){ //TODO
        eff.SetCentBinOut("0 100");
    }else{
        eff.SetCentBinOut("0 10 20 40 60 90");
    }
    //eff.SetCentBinOut("0 100"); //TODO for pA
    eff.fRT = &rt;
    eff.LoadHist();
    eff.EvaluateAll();
    eff.Write();
}
