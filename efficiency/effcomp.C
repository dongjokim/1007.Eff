int qColor[] = { kBlack, kRed, kBlue, kMagenta, kGreen, kViolet, kOrange };
int qMarker[]= { 29,30,20, 24, 21, 25, 33, 27, 34, 28 };
int qLineStyle[] = { 1, 2, 5, 8, 9};

enum { kJFilterTPCOnly, kJFilterRaa, kJFilterTightDCA, kJFilterGlobalDCA, kJFilterGlobalSDD, kJFilterHybrid,
    kJNFilter
};
TString strTrk[] = { "TPCOnly","Raa","Global+TightDCA","Global+DCA", "Global+SDD", "Hybrid",
};

TString labelStr[] = {"",""};

void effcomp(TString ref="eff/EffLHC17pFast-LHC17l3b-0-LHC17l3b.root",TString infile="eff/Effcent_woSDD-LHC17l3b-0-LHC17l3b.root") {

	int iC = 0;
	const int Nsets = 2;
	TFile *fin[Nsets];
	fin[0] = new TFile(ref.Data(),"open");
	fin[1] = new TFile(infile.Data(),"open");
	TGraphErrors * gCor[Nsets][kJNFilter] = {0x0};

	for(int i=0;i<Nsets;i++) {
		for(int icut=0;icut<kJNFilter;icut++) {
			gCor[i][icut] = (TGraphErrors*)fin[i]->Get(Form("Efficiency/gCor00%02d%02d",iC,icut));
			//if(i==0) gCor[i][icut] = (TGraphErrors*)fin[i]->Get(Form("Efficiency/gCor00%02d%02d",iC,icut));
			//if(i==1) gCor[i][icut] = (TGraphErrors*)fin[i]->Get(Form("Eff/gCor00%02d%02d",iC,icut));
		}
	}

	cout <<"Drawing "<< endl;
    mc(1,1.1);
    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(0);mpad->SetGridx(0);
    mpad->SetLogy(0);mpad->SetGridy(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 0, 20, 10, 0.3, 1.15);
    hset( *hfr, "p_{T}","Correction Factor");
    hfr->Draw();

    TString name = Form("pp: 5TeV %s", strTrk[kJFilterGlobalSDD].Data());

    leg = new TLegend(0.25,0.2,0.6,0.4,name.Data(),"brNDC");
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.04);
    leg->AddEntry((TObject*)0, name, "");
    for(int i=0;i<Nsets;i++) {
	    gCor[i][kJFilterGlobalSDD]->SetMarkerStyle(qMarker[i]);
	    gCor[i][kJFilterGlobalSDD]->SetMarkerColor(qColor[i]);
	    gCor[i][kJFilterGlobalSDD]->Draw("p,same");
	    //leg->AddEntry(gCor[i][kJFilterGlobalSDD], Form("%s %s", strTrk[kJFilterGlobalSDD].Data(), fin[i]->GetName()), "p");
	    leg->AddEntry(gCor[i][kJFilterGlobalSDD], Form("%s %s", strTrk[kJFilterGlobalSDD].Data(), labelStr[i].Data()), "p");
    }
    leg->Draw();

    TGraphErrors * ratio = get_ratio( gCor[0][kJFilterGlobalSDD], gCor[1][kJFilterGlobalSDD] );

    printGrr(gCor[0][kJFilterGlobalSDD]);




    mc(2,1.1);
    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(0);mpad->SetGridx(0);
    mpad->SetLogy(0);mpad->SetGridy(0);
    gPad->SetLeftMargin(0.17);
    hfr = new TH2F("hfr"," ", 10, 0, 20, 10, 0.93, 1.07);
    hset( *hfr, "p_{T}","Ratio");
    hfr->Draw();

    leg = new TLegend(0.2,0.2,0.7,0.5,"","brNDC");
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.04);
    leg->AddEntry((TObject*)0, name, "");
    ratio->SetMarkerStyle(qMarker[0]);
    ratio->SetMarkerColor(qColor[0]);
    ratio->Draw("p,same");
    leg->Draw();


}

TGraphErrors* get_ratio( TGraphErrors * l, TGraphErrors *r ){
	TGraphErrors * gr_ratio = new TGraphErrors( l->GetN() );
	TGraph ger( r->GetN(),  r->GetX(), r->GetEY() );
	for( int i=0; i< l->GetN(); i++ ){
		double x = l->GetX()[i];
		double y1 = l->GetY()[i];
		double ey1 = l->GetEY()[i];
		double y2 = r->Eval(x);
		double ey2 = ger.Eval(x);


		double ratio = y1 / y2;
		gr_ratio->SetPoint( i,  x, ratio);
		gr_ratio->SetPointError( i,  0, ratio*TMath::Sqrt( ey1*ey1/y1/y1+ey2*ey2/y2/y2));
	}
	return gr_ratio;
}
