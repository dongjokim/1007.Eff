int qColor[] = { kBlack, kRed, kBlue, kMagenta, kGreen, kViolet, kOrange, kPink};
//int qMarker[]= { 29,30,20, 24, 21, 25, 33, 27, 34, 28 };
int qMarker[]= { 24,24,24, 24, 21, 25, 33, 27, 34, 28 };
int qLineStyle[] = { 1, 2, 5, 8, 9};
TH2F * hfr;
TLegend *led;
int icolor;
int iCan;
TString Name;
TString FigDir;
TString OutTag;
void DrawEffTrigg(TString input, int DrawMode = 1 ){

    TFile * inRoot = new TFile(input, "READ");
    if( !inRoot ){
        cout<<"NO Files : "<<input<<endl;
        gApplication->Exit(1);
    }

    FigDir = "figs/";
    OutTag = input;
    OutTag.ReplaceAll("/","_");
    OutTag.ReplaceAll(".","");
    OutTag.ReplaceAll("root","");
    OutTag = FigDir+"/"+OutTag;

    TDirectory *effD = (TDirectory*) inRoot->Get("Efficiency");
    effD->cd();

    vector<double> CentBin = GetVecD("CentralityBin", 1);
    cout << "---------------------" << endl;
    cout << CentBin.size()-1 << endl;
    vector<double> VtxBin = GetVecD("VtxBin",1);
    vector<double> PtBin = GetVecD("PtBin",1);
    vector<TString> TrackCut = GetVecS("TrackCut",1);
    /* 
    TNamed *Period = (TNamed*) effD->Get("Period");
    if(Period)Period->Print();
    TNamed *Runnumber = (TNamed*) effD->Get("RunNumber");
    if(Runnumber) Runnumber->Print();
    */

    if( VtxBin.size() == 0 ){
        VtxBin.push_back(-10);
        VtxBin.push_back(10);
        TrackCut.push_back("TPCOnly"); 
        TrackCut.push_back("Raa"); 
        TrackCut.push_back("Global+Tight DCA"); 
        TrackCut.push_back("Global+DCA"); 
        TrackCut.push_back("Global+SDD"); 
        TrackCut.push_back("Hybrid");
    }


    gROOT->LoadMacro("AliJRunTable.cxx+g");
    AliJRunTable rt;
    //rt.SetPeriod( Period->GetTitle() );
    rt.ParseString( input ); // For Back compat
    TString Title = Form("%s %s %s",rt.GetPeriodName().Data(), rt.GetBeamTypeStr().Data(), rt.GetEnergyStr());
    cout<<"=== Title "<<Title<<endl;



    const int kNVtx = 2;
    const int kNCent = 20;
    const int kNCut = 10;
    TGraphErrors *gCor[kNVtx][kNCent][kNCut] = {0};
    TGraphErrors *gEff[kNVtx][kNCent][kNCut] = {0};
    TGraphErrors *gCon[kNVtx][kNCent][kNCut] = {0};
    TGraphErrors *gEffTrigg[kNVtx][kNCent][kNCut] = {0};
    TGraphErrors *gEffTriggVtx[kNVtx][kNCent][kNCut] = {0};
    TString CentStr[kNCent];

    //==== READ Graph
    for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
        for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
            if( rt.IsPP() ){
                CentStr[icent] = " ";
            }else{
                CentStr[icent] = Form(" Cent %.0f~%.0f",CentBin[icent], CentBin[icent+1]);
            }
            for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                gCor[ivtx][icent][icut] = (TGraphErrors*)effD->Get(Form("gCor%02d%02d%02d",ivtx,icent,icut));
                gEff[ivtx][icent][icut] = (TGraphErrors*)effD->Get(Form("gEff%02d%02d%02d",ivtx,icent,icut));
                gCon[ivtx][icent][icut] = (TGraphErrors*)effD->Get(Form("gCon%02d%02d%02d",ivtx,icent,icut));
                gEffTrigg[ivtx][icent][icut] = (TGraphErrors*)effD->Get(Form("gEffTrigg%02d%02d%02d",ivtx,icent,icut));
                gEffTriggVtx[ivtx][icent][icut] = (TGraphErrors*)effD->Get(Form("gEffTriggVtx%02d%02d%02d",ivtx,icent,icut));

            }
        }
    }
    TGraphErrors *ge = 0;
    iCan = 0;

    double lpt = 1e-1;
    double hpt = 100;
    double ly = 0.0;
    double hy = 1.1;

    if( DrawMode & 1 ){
        //==== DRAW Cor
        for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
            for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
                mc(++iCan);
                hset0( Form("Cor_cent%d",icent),lpt,hpt,ly,hy, "p_{T} [GeV/c]", "Correction Factor");
                leg0(0.7,0.6,0.9,0.9,Title);
                if(!rt.IsPP()) led->AddEntry(hfr, CentStr[icent], "p" );
                for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                    //if( icut!=0 && icut!=1 && icut!=3 && icut!=7 ) continue;
                    ge = gCor[ivtx][icent][icut];
                    DrawGraph( ge, icolor++, TrackCut[icut] );
                }
                FinishMC();

            }
        }

        //==== DRAW Eff 
        for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
            for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
                mc(++iCan);
                hset0( Form("Eff_cent%d",icent),lpt,hpt,ly,hy, "p_{T} [GeV/c]", "Pure Efficiency");
                leg0(0.7,0.6,0.9,0.9,Title);
                if(!rt.IsPP()) led->AddEntry(hfr, CentStr[icent], "p" );
                for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                    //if( icut!=0 && icut!=1 && icut!=3 && icut!=7 ) continue;
                    ge = gEff[ivtx][icent][icut];
                    DrawGraph( ge, icolor++, TrackCut[icut] );
                }
                FinishMC();
            }
        }

        //==== DRAW Con
        for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
            for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
                mc(++iCan);
                hset0( Form("Con_cent%d",icent),lpt,hpt,0,1, "p_{T} [GeV/c]", "Contamination");
                leg0(0.7,0.6,0.9,0.9,Title);
                if(!rt.IsPP()) led->AddEntry(hfr, CentStr[icent], "p" );
                for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                    //if( icut!=0 && icut!=1 && icut!=3 && icut!=7 ) continue;
                    ge = gCon[ivtx][icent][icut];
                    DrawGraph( ge, icolor++, TrackCut[icut] );
                }
                FinishMC();
            }
        }

        //==== DRAW EffTrigg
        for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
            for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
                mc(++iCan);
                hset0( Form("EffTrigg_cent%d",icent),lpt,hpt,0.9,1.1, "p_{T} [GeV/c]", "Trigg eff of tracks");
                leg0(0.7,0.6,0.9,0.9,Title);
                if(!rt.IsPP()) led->AddEntry(hfr, CentStr[icent], "p" );
                for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                    //if( icut!=0 && icut!=1 && icut!=3 && icut!=7 ) continue;
                    ge = gEffTrigg[ivtx][icent][icut];
                    DrawGraph( ge, icolor++, TrackCut[icut] );
                }
                FinishMC();
            }
        }
        //==== DRAW EffTriggVtx
        for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
            for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
                mc(++iCan);
                hset0( Form("EffTriggVtx_cent%d",icent),lpt,hpt,0.1,1.1, "p_{T} [GeV/c]", "TriggVtx eff of tracks");
                leg0(0.7,0.6,0.9,0.9,Title);
                if(!rt.IsPP()) led->AddEntry(hfr, CentStr[icent], "p" );
                for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                    //if( icut!=0 && icut!=1 && icut!=3 && icut!=7 ) continue;
                    ge = gEffTriggVtx[ivtx][icent][icut];
                    DrawGraph( ge, icolor++, TrackCut[icut] );
                }
                FinishMC();
            }
        }
    }


    if( DrawMode & 2 ){
        //==== DRAW Cor
        for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
            for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                //if( icut!=0 && icut!=1 && icut!=3 && icut!=7 ) continue;
                mc(++iCan);
                hset0( Form("Cor_cut%d",icut),lpt,hpt,ly,hy, "p_{T} [GeV/c]", "Correction Factor");
                leg0(0.7,0.6,0.9,0.9,Title);
                led->AddEntry(hfr, TrackCut[icut] , "");
                for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
                    ge = gCor[ivtx][icent][icut];
                    DrawGraph( ge, icolor++, CentStr[icent] );
                }
                FinishMC();
            }
        }

        //==== DRAW Eff 
        for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
            for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                //if( icut!=0 && icut!=1 && icut!=3 && icut!=7 ) continue;
                mc(++iCan);
                hset0( Form("Eff_Cut%d",icut),lpt,hpt,ly,hy, "p_{T} [GeV/c]", "Pure Efficiency");
                leg0(0.7,0.6,0.9,0.9,Title);
                led->AddEntry(hfr, TrackCut[icut] , "");
                for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
                    ge = gEff[ivtx][icent][icut];
                    DrawGraph( ge, icolor++, CentStr[icent] );
                }
                FinishMC();
            }
        }

        //==== DRAW Con
        for( UInt_t ivtx=0;ivtx<VtxBin.size()-1;ivtx++ ){
            for( UInt_t icut=0;icut<TrackCut.size();icut++ ){
                //if( icut!=0 && icut!=1 && icut!=3 && icut!=7 ) continue;
                mc(++iCan);
                hset0( Form("Con_Cut%d",icut),lpt,hpt,0,0.5, "p_{T} [GeV/c]", "Contamination");
                leg0(0.7,0.6,0.9,0.9,Title);
                led->AddEntry(hfr, TrackCut[icut] , "");
                for( UInt_t icent=0;icent<CentBin.size()-1;icent++ ){
                    ge = gCon[ivtx][icent][icut];
                    DrawGraph( ge, icolor++, CentStr[icent] );
                }
                FinishMC();
            }
        }
    }

}



void hset0(TString name, double x0, double x1, double y0, double y1, TString xtitle, TString ytitle){
    Name=name;
    hfr = new TH2F(Form("hfr%d",iCan), "hfr", 10, x0, x1, 10, y0, y1);
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    gPad->SetLeftMargin(0.17);
    gPad->SetGridx(1);
    gPad->SetGridy(2);
    gPad->SetLogx(1);
    hset(*hfr, xtitle, ytitle,0.9,1.4, 0.05,0.05, 0.01,0.001, 0.04,0.03, 505,510 );
    hfr->Draw();
    icolor=0;
}
void leg0(double x0, double x1, double y0, double y1, TString legstr ){
    led = new TLegend(x0, x1, y0, y1, legstr.Data());
    led->SetTextSize(0.04);
}

void FinishMC(){
    led->Draw();

    gPad->GetCanvas()->SaveAs( TString(OutTag+"_"+Name+".pdf").Data() );
}

void DrawGraph( TGraphErrors * gr, int i, TString legstr ){
    if( !gr ) return;
    gr->SetLineColor(qColor[i]);
    gr->SetMarkerStyle(qMarker[i]);
    gr->SetMarkerColor(qColor[i]);
    gr = (TGraphErrors*)gr->DrawClone("samep");
    led->AddEntry(gr, legstr, "p");
}

vector<double> GetVecD(TString name, UInt_t mode ){
    TAxis *axis = (TAxis*)gDirectory->Get(name);
    vector<double> v;
    if( ! axis ) {
        cout << "DEBUG: sizeof vec=" << v.size() << endl;
        return v;
    }
    for( UInt_t i=1;i<=axis->GetLast()+1;i++ ){
        v.push_back(axis->GetBinLowEdge(i));
    }
    if( mode ){
        cout<<"==== "<<name<<" : "<<v.size()<<endl;
        for( UInt_t i=0;i<v.size();i++){
            cout<<v[i]<<", ";
        }
        cout<<endl;
    }
    return v;
}

vector<TString> GetVecS(TString name, UInt_t mode ){
    TAxis *axis = (TAxis*)gDirectory->Get(name);
    vector<TString> v;
    if( ! axis ) return v;
    for( UInt_t i=1;i<=axis->GetLast();i++ ){
        v.push_back(axis->GetBinLabel(i));
    }
    if( mode ){
        cout<<"==== "<<name<<" : "<<v.size()<<endl;
        for( UInt_t i=0;i<v.size();i++){
            cout<<v[i]<<" ";
        }
        cout<<endl;
    }
    return v;
}

