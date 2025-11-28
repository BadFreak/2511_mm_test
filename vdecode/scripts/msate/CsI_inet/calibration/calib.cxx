std::vector<float> findPeak(TH1D* h){
	std::vector<float> ret;
	for(int i=0;i<h->GetNbinsX();i++){
		if(h->GetBinContent(i+1)>100.){
			float mean = h->GetBinCenter(i+1);
			TF1 *f1=new TF1("f1","[0]*TMath::Gaus(x,[1],[2])",mean-200,mean+200);
			f1->SetParameters(100,mean,100);
			f1->SetParLimits(2,30,100);
			h->Fit("f1","+q","",mean-200,mean+200);
			float output_mean = f1->GetParameter(1);
			float output_rms = f1->GetParameter(2);
			float output_mag = f1->GetParameter(0);
			std::cout<<output_mean<<" "<<output_rms<<" "<<output_mag<<std::endl;
			ret.emplace_back(output_mean);
			while(h->GetBinContent(i+1)>50.)i++;
			delete f1;
		}
	}
	return ret;
}
void calib(){
	TGraph *gh[4],*gl[4];
	auto f = TFile::Open("hist.root","READ");
	int bhchns[4]={2,3,4,9};
	int blchns[4]={4,9,14,19};
	for(int i=0;i<4;i++){
		gh[i] = new TGraph();
		gh[i]->SetName(Form("gh_%d",bhchns[i]));
		gl[i] = new TGraph();
		gl[i]->SetName(Form("gl_%d",blchns[i]));
		auto hbh = (TH1D*)f->Get(Form("hbh_%d",bhchns[i]));
		auto hbl = (TH1D*)f->Get(Form("hbl_%d",blchns[i]));
		auto hpeaks = findPeak(hbh);
		auto lpeaks = findPeak(hbl);
		int index=0;
		for(auto h:hpeaks){
			gh[i]->SetPoint(index,index+1,h);
			index++;
		}
		index=0;
		for(auto l:lpeaks){
			gl[i]->SetPoint(index,index+1,l);
			index++;
		}
	}
	std::cout<<"Peaks found"<<std::endl;
	auto fout = new TFile("peaks.root","RECREATE");
	fout->cd();
	for(int i=0;i<4;i++){
		gh[i]->Write();
		gl[i]->Write();
	}
	fout->Close();
}
