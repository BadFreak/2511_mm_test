#include "ComReader.hh"

bool ComReader::openFile(const std::string& filename){
    file = new std::ifstream(filename,std::ios::binary);
    if (!file->is_open()) {
        std::cerr << "Can not open file: "<<filename << std::endl;
        return false;
    }
    return true;
}

bool ComReader::findHead()
{
	char buffer[1];
	while (file->rdbuf()->sgetn(buffer, sizeof(buffer)) > 0)
	{
		unsigned char byte = static_cast<unsigned char>(buffer[0]);
		if (byte == 0x55)
		{
			char buffer2[3];
			if (file->rdbuf()->sgetn(buffer2, sizeof(buffer2)) > 0)
			{
				unsigned char b1 = static_cast<unsigned char>(buffer2[0]);
				unsigned char b2 = static_cast<unsigned char>(buffer2[1]);
				unsigned char b3 = static_cast<unsigned char>(buffer2[2]);
				if (b1 == 0xaa && b2 == 0xeb && b3 == 0x90)
				{
					char buffer_length[4];
					file->rdbuf()->sgetn(buffer_length, sizeof(buffer_length));
					unsigned char length_0 = static_cast<unsigned char>(buffer_length[2]);
					unsigned char length_1 = static_cast<unsigned char>(buffer_length[3]);
					PackageLength = static_cast<unsigned short>(length_0) << 8 | length_1;
					// std::cout << "PackageLength: " << PackageLength << std::endl;
					if (PackageLength == 116 || PackageLength == 44 || PackageLength == 6668 || PackageLength == 2060)
					{ // Calo or CsITK
						file->seekg(-4, std::ios::cur);
						nhead++;
						if (!readFEE())
						{
							logger->error("Failed to read FEE data");
						}
					}
					else
					{
						// // Judge if it is HK
						// unsigned char b2 = static_cast<unsigned char>(buffer_length[0]);
						// if ((b2 & 0x0F) == 0x0B)
						// { // Temperature data packages
						// 	file->seekg(-3, std::ios::cur);
						// 	unsigned char x = (unsigned char)((b2 >> 4) & 0x0F);
						// 	FEEID = int(x) - 4;
						// 	if (FEEID >= 1 && FEEID <= 4)
						// 	{
						// 		nHKHead++;
						// 		if (!readHK())
						// 		{
						// 			logger->error("Failed to read HK data");
						// 		}
						// 	}
						// 	else{
						// 		if(FEEID!=5 && FEEID!=-4 && FEEID!=-3)std::cout << "Invalid FEEID for HK: " << FEEID << std::endl;
						// 	}
						// }
						// else
						// {
						// 	// std::cout << "Not HK Package with length : "<<PackageLength << std::endl;
						// }
					}
				}
				else
				{
					file->seekg(-1, std::ios::cur); // Not a scientific data head, seek back
				}
			}
			else
			{
				file->seekg(-3, std::ios::cur);
				return false;
			}
		}
	}
	return false;
}

bool ComReader::readFEE(){ // Read length, FEEID and other info
	char buffer[12];
	file->rdbuf()->sgetn(buffer,sizeof(buffer));
	auto length_0 = static_cast<unsigned char>(buffer[2]);
	auto length_1 = static_cast<unsigned char>(buffer[3]);
	PackageLength = static_cast<unsigned short>(length_0)<<8 | length_1;
	// std::cout << "PackageLength: " << PackageLength << std::endl;
	unsigned char bfee = static_cast<unsigned char>(buffer[11]);
	PackageID = static_cast<int>(bfee & 0x0F); //Package class
	FEEID = static_cast<int>((bfee & 0xF0) >> 4); //FEEID 前四位是FEEID 后四位是包类型标识
	if(FEEID >4 && FEEID < 9){ // Calo
		nCaloHead++;
		if(PackageLength == 6668){ // Calo wave
			if(!readCaloWave()){
				logger->error("Failed to read Calo wave data returning to findHead");
			}
		}
		else if(PackageLength == 116){
			if(!readCalo()){
				logger->error("Failed to read Calo data");
			}
		}
		else{
			logger->error("Unexpected Calo package length: {}", PackageLength);
			return false;
		}
	}
	else if(FEEID == 2){ // CsI
		nCsIHead++;
		if(PackageLength == 2060){ // CsI wave
			if(!readCsIWave()){
				logger->error("Failed to read CsI wave data");
			}
		}
		else if(PackageLength == 44){
			if(!readCsI()){
				logger->error("Failed to read CsI data");
			}
		}
		else{
			logger->error("Unexpected CsI package length: {}", PackageLength);
			return false;
		}
	}
	else{
		// logger->info("Unknown FEEID: {}", FEEID);
	}
	return true;
}

void ComReader::clear(){
	CellID.clear();
	CellADC.clear();
	CellPLAT.clear();
}

void ComReader::clearHK(){
	C0.clear();
	C1.clear();
	C2.clear();
	T0.clear();
	T1.clear();
	T2.clear();
	T3.clear();
	C0=std::vector<float>(4,0.);
	C1=std::vector<float>(4,0.);
	C2=std::vector<float>(4,0.);
	T0=std::vector<float>(4,0.);
	T1=std::vector<float>(4,0.);
	T2=std::vector<float>(4,0.);
	T3=std::vector<float>(4,0.);
}

bool ComReader::readCsI(){
	if(PackageLength!=44){
		logger->error("Unexpected CsI package length: {}", PackageLength);
		return false;
	}
	char buf[38];
	file->rdbuf()->sgetn(buf,sizeof(buf));
	//Trigger ID
	unsigned char triggerid_0 = static_cast<unsigned char>(buf[32]);
	unsigned char triggerid_1 = static_cast<unsigned char>(buf[33]);
	TriggerID_csi = static_cast<unsigned short>(triggerid_0 & 0x0F)<<8 | triggerid_1;
	TriggerIDMM_csi = readAmpTriggerIDMM();
	readBuffer(buf, 8);
	csiTree->Fill();
	clear();
	CsI_EventID++;
	return true;
}

bool ComReader::readCalo(){
	if(PackageLength!=116){
		logger->error("Unexpected Calo package length: {}", PackageLength);
		return false;
	}
	char buf[110];
	file->rdbuf()->sgetn(buf,sizeof(buf));
	//Trigger ID
	unsigned char triggerid_0 = static_cast<unsigned char>(buf[104]);
	unsigned char triggerid_1 = static_cast<unsigned char>(buf[105]);
	TriggerID = static_cast<unsigned short>(triggerid_0 & 0x0F)<<8 | triggerid_1;
	readBuffer(buf, 26);
	Calo_EventCount++;
	TriggerIDMM = readAmpTriggerIDMM();
	if(Calo_EventCount==4){
		caloTree->Fill();
		clear();
		Calo_EventID++;
		Calo_EventCount=0;
	}
	return true;
}

int ComReader::readAmpTriggerIDMM(){
	char buf[6];
	file->rdbuf()->sgetn(buf,sizeof(buf));
	// triggerIDMM consist of last 4 bytes
	unsigned char triggerid_0 = static_cast<unsigned char>(buf[2]);
	unsigned char triggerid_1 = static_cast<unsigned char>(buf[3]);
	unsigned char triggerid_2 = static_cast<unsigned char>(buf[4]);
	unsigned char triggerid_3 = static_cast<unsigned char>(buf[5]);
	TriggerIDMM = static_cast<unsigned int>(triggerid_0) << 24 | static_cast<unsigned int>(triggerid_1) << 16 | static_cast<unsigned int>(triggerid_2) << 8 | static_cast<unsigned int>(triggerid_3);
	std::cout << std::hex << "MM TriggerID: " << TriggerIDMM << std::endl;
	return TriggerIDMM;
} 

int ComReader::readWaveTriggerIDMM(){
	char buf[4];
	file->rdbuf()->sgetn(buf,sizeof(buf));
	// triggerIDMM consist of 4 byte
	unsigned char triggerid_0 = static_cast<unsigned char>(buf[0]);
	unsigned char triggerid_1 = static_cast<unsigned char>(buf[1]);
	unsigned char triggerid_2 = static_cast<unsigned char>(buf[2]);
	unsigned char triggerid_3 = static_cast<unsigned char>(buf[3]);
	TriggerIDMM = static_cast<unsigned int>(triggerid_0) << 24 | static_cast<unsigned int>(triggerid_1) << 16 | static_cast<unsigned int>(triggerid_2) << 8 | static_cast<unsigned int>(triggerid_3);
	// std::cout << std::dec << "MM TriggerID: " << TriggerIDMM << std::endl;
	return TriggerIDMM;
}

bool ComReader::readCaloWave(){
	// logger->info("Reading Calo wave data");
	if(PackageLength!=6668){
		logger->error("Unexpected Calo wave package length: {}", PackageLength);
		return false;
	}
	char buf[6664];
	file->rdbuf()->sgetn(buf,sizeof(buf));
	//Last two bytes should be 5aa5
	unsigned char b5a = static_cast<unsigned char>(buf[6662]);
	unsigned char ba5 = static_cast<unsigned char>(buf[6663]);
	if(b5a != 0x5a || ba5 != 0xa5){
		logger->error("Unexpected Calo wave package end: {}", b5a << 8 | ba5);
		return false;
	}
	// 26 channels, 128 sample points 
	for(int chn_i = 0; chn_i < 26; chn_i++){
		std::vector<unsigned short> channel_buffer;
		for(int sample_i = 0; sample_i < 128; sample_i++){
			auto d1 = static_cast<unsigned char>(buf[chn_i * 128 * 2 + sample_i * 2]);
			auto d2 = static_cast<unsigned char>(buf[chn_i * 128 * 2 + sample_i * 2 + 1]);
			unsigned short value = (static_cast<unsigned short>(d1) << 8) | d2;
			channel_buffer.emplace_back(value);
		}
		//Pedestal is the mean of the first 16 samples
		int pedestal = std::accumulate(channel_buffer.begin(), channel_buffer.begin() + 16, 0.) / 16.;
		//Amplitude is the maximum value of the channel
		int amplitude = *std::max_element(channel_buffer.begin(), channel_buffer.end());
		int feeid=FEEID-4;
		std::pair<int, int> gid_cryid = std::pair<int, int>(0, 0);
		if (channelMap_Calo.count(std::pair<int, int>(feeid, chn_i)) == 0)
		{
			if (feeid == 3 || feeid == 4)
			{
				if (chn_i == 12 || chn_i == 25)
				{
					continue;
				}
			}
		}
		gid_cryid = channelMap_Calo[std::pair<int, int>(feeid, chn_i)];
		int tmp_cellid = gid_cryid.second * 100000 + 10000 * feeid + 1000 * (feeid % 2) + 100 * (gid_cryid.first) + chn_i;
		CellID.emplace_back(tmp_cellid);
		CellADC.emplace_back(amplitude);
		CellPLAT.emplace_back(pedestal);
	}
	//Trigger ID
	unsigned char triggerid_0 = static_cast<unsigned char>(buf[6656]);
	unsigned char triggerid_1 = static_cast<unsigned char>(buf[6657]);
	TriggerID = static_cast<unsigned short>(triggerid_0 & 0x0F)<<8 | triggerid_1;
	Calo_EventCount++;
	//Trigger ID MicroMegas 
	TriggerIDMM = readWaveTriggerIDMM();
	std::cout << "TriggerIDMM: " << TriggerIDMM << std::endl;
	if(Calo_EventCount==4){
		caloTree->Fill();
		clear();
		Calo_EventID++;
		Calo_EventCount=0;
	}
	return true;
}

bool ComReader::readCsIWave(){
	// logger->info("Reading CsI wave data");
	if(PackageLength!=2060){
		logger->error("Unexpected CsI wave package length: {}", PackageLength);
		return false;
	}
	char buf[2056];
	file->rdbuf()->sgetn(buf,sizeof(buf));
	//Last two bytes should be 5aa5
	unsigned char b5a = static_cast<unsigned char>(buf[2054]);
	unsigned char ba5 = static_cast<unsigned char>(buf[2055]);
	if(b5a != 0x5a || ba5 != 0xa5){
		logger->error("Unexpected CsI wave package end: {}", b5a << 8 | ba5);
		return false;
	}
	// 8 CsI channels, 128 sample points 
	for(int chn_i = 0; chn_i < 8; chn_i++){
		std::vector<unsigned short> channel_buffer;
		for(int sample_i = 0; sample_i < 128; sample_i++){
			auto d1 = static_cast<unsigned char>(buf[chn_i * 128 * 2 + sample_i * 2]);
			auto d2 = static_cast<unsigned char>(buf[chn_i * 128 * 2 + sample_i * 2 + 1]);
			unsigned short value = (static_cast<unsigned short>(d1) << 8) | d2;
			channel_buffer.emplace_back(value);
		}
		//Pedestal is the mean of the first 16 samples
		int pedestal = std::accumulate(channel_buffer.begin(), channel_buffer.begin() + 16, 0.) / 16.;
		//Amplitude is the maximum value of the channel
		int amplitude = *std::max_element(channel_buffer.begin(), channel_buffer.end());
		int feeid=FEEID-4;
		int tmp_cellid = (chn_i+1)*100000;
		CellID.emplace_back(tmp_cellid);
		CellADC.emplace_back(amplitude);
		CellPLAT.emplace_back(pedestal);
	}
	//Trigger ID	
	unsigned char triggerid_0 = static_cast<unsigned char>(buf[2048]);
	unsigned char triggerid_1 = static_cast<unsigned char>(buf[2049]);
	TriggerID = static_cast<unsigned short>(triggerid_0 & 0x0F)<<8 | triggerid_1;
	TriggerIDMM_csi = readWaveTriggerIDMM();
	csiTree->Fill();
	clear();
	CsI_EventID++;
	return true;
}

bool ComReader::readBuffer(char *b,const int& NCHN) {
	if(NCHN==26){ //Calo
		for(size_t chn_i = 0; chn_i < NCHN; chn_i++)
		{
			auto d1 = static_cast<unsigned char>(b[chn_i * 4]);
			auto d2 = static_cast<unsigned char>(b[chn_i * 4 + 1]);
			auto d3 = static_cast<unsigned char>(b[chn_i * 4 + 2]);
			auto d4 = static_cast<unsigned char>(b[chn_i * 4 + 3]);
			unsigned short plat = (static_cast<unsigned short>(d1) << 8) | d2;
			unsigned short maxi = (static_cast<unsigned short>(d3) << 8) | d4;
			std::pair<int, int> gid_cryid = std::pair<int, int>(0, 0);
			int feeid=FEEID-4;
			if (channelMap_Calo.count(std::pair<int, int>(feeid, chn_i)) == 0)
			{
				if (feeid == 3 || feeid == 4)
				{
					if (chn_i == 12 || chn_i == 25)
					{
						continue;
					}
				}
			}
			gid_cryid = channelMap_Calo[std::pair<int, int>(feeid, chn_i)];
			if (gid_cryid.second == 0)
			{
				std::cerr << "No channel map for FEEID: " << feeid << " chn: " << chn_i << std::endl;
			}
			int tmp_cellid = gid_cryid.second * 100000 + 10000 * feeid + 1000 * (feeid % 2) + 100 * (gid_cryid.first) + chn_i;
			CellID.emplace_back(tmp_cellid);
			CellADC.emplace_back(static_cast<int>(maxi));
			CellPLAT.emplace_back(plat);
		}
	}
	else if(NCHN==8){//CsI
		for(int chn_i = 0; chn_i < NCHN; chn_i++)
		{
			auto d1 = static_cast<unsigned char>(b[chn_i * 4]);
			auto d2 = static_cast<unsigned char>(b[chn_i * 4 + 1]);
			auto d3 = static_cast<unsigned char>(b[chn_i * 4 + 2]);
			auto d4 = static_cast<unsigned char>(b[chn_i * 4 + 3]);
			unsigned short plat = (static_cast<unsigned short>(d1) << 8) | d2;
			unsigned short maxi = (static_cast<unsigned short>(d3) << 8) | d4;
			int tmp_cellid = (chn_i+1)*100000;
			CellID.emplace_back(tmp_cellid);
			CellADC.emplace_back(static_cast<int>(maxi));
			CellPLAT.emplace_back(plat);
		}
	}
	else{
		std::cerr << "Unknown channel number: " << NCHN << std::endl;
	}
	return true;
}

bool ComReader::readHK(){
	if(HK_EventCount==0){
		TPointID++;
		clearHK();
	}
	//Start to read
	char rest[60];
	file->rdbuf()->sgetn(rest,sizeof(rest));
	auto f_current = [this,rest](const int& i){
		int index = i*2;
		float cif=-1.;
		auto ci0 = static_cast<unsigned char>(rest[index]);
		auto ci1 = static_cast<unsigned char>(rest[index+1]);
		unsigned short ci = (static_cast<unsigned short>(ci0 & 0x0F) << 8) | ci1;
		if(i==1){
			cif = (ci*2500./4096.-245.1)/20./0.01;
			// if(dettype == DetType::ITK) cif/=4.;
		}
		else{
			cif = (ci*2500./4096.-245.1)/20./0.125;
		}
		return cif;
	};
	C0[FEEID-1] = f_current(1);
	C1[FEEID-1] = f_current(2);
	C2[FEEID-1] = f_current(3);
	// if(dettype == DetType::ITK && FEEID==2){
	// 	C1[1]/=2.;
	// 	C2[1]/=2.;
	// }
	// else if(dettype == DetType::CALO){
	// 	if(FEEID==3)C1[2]/=2.;
	// 	if(FEEID==4)C1[3]/=2.;
	// }
	//Temperature
	auto f_temperature = [rest](const int& i){
		int index = 6 + i*2;
		float tif=-1.;
		auto ti0 = static_cast<unsigned char>(rest[index]);
		auto ti1 = static_cast<unsigned char>(rest[index+1]);
		unsigned short ti = (static_cast<unsigned short>(ti0 & 0x0F) << 8) | ti1;
		float fx = 10.*(ti*2500./4096.*2000./2100.)/(2500.-(ti*2500./4096.*2000./2100.));
		tif = (17.-fx)/(12./25.);
		return tif;
	};
	T0[FEEID-1] = f_temperature(1);
	T1[FEEID-1] = f_temperature(2);
	T2[FEEID-1] = f_temperature(3);
	T3[FEEID-1] = f_temperature(4);
	HK_EventCount++;
	if(HK_EventCount==4){
		hkTree->Fill();
		HK_EventCount=0;
	}
	return true;
}

void ComReader::initCaloTree(){ //Initialize Calo Tree
	caloTree = new TTree("caloTree","caloTree");
	caloTree->Branch("PackageID",&PackageID);
	caloTree->Branch("TriggerID",&TriggerID);
	caloTree->Branch("EventID",&Calo_EventID);
	caloTree->Branch("CellID",&CellID);
	caloTree->Branch("CellADC",&CellADC);
	caloTree->Branch("CellPLAT",&CellPLAT);
	// TriggerIDMM 
	caloTree->Branch("TriggerIDMM",&TriggerIDMM);
}

void ComReader::initCsITree(){ //Initialize CsITK Tree
	csiTree = new TTree("csiTree","csiTree");
	csiTree->Branch("PackageID",&PackageID);
	csiTree->Branch("TriggerID",&TriggerID_csi);
	csiTree->Branch("EventID",&CsI_EventID);
	csiTree->Branch("CellID",&CellID);
	csiTree->Branch("CellADC",&CellADC);
	csiTree->Branch("CellPLAT",&CellPLAT);
	// TriggerIDMM_csi
	csiTree->Branch("TriggerIDMM",&TriggerIDMM_csi);
}

void ComReader::initHKTree(){ //Initialize HK Tree
	hkTree = new TTree("hkTree","hkTree");
	hkTree->Branch("TPoint",&TPointID);
	hkTree->Branch("C0",&C0);
	hkTree->Branch("C1",&C1);
	hkTree->Branch("C2",&C2);
	hkTree->Branch("T0",&T0);
	hkTree->Branch("T1",&T1);
	hkTree->Branch("T2",&T2);
	hkTree->Branch("T3",&T3);
}

void ComReader::decode(const std::string& filename){
	logger->info("Decoding file: {}", filename);
	getOutputName(std::string("result_"),filename);
    fout = new TFile(oname.c_str(),"RECREATE");
	logger->info("Initialized output file: {}", oname);
	initCaloTree();
	initCsITree();
	initHKTree();
	logger->info("Initialized output trees");
    openFile(filename);
	auto total_size = file->seekg(0, std::ios::end).tellg();
	logger->info("Total file size: {}", int(total_size));
	file->seekg(0, std::ios::beg);
    while(true){
        if(!findHead()){
			break;
		}
    }

	logger->info("Total head: {}", nhead);
	logger->info("CALO: {} expected {}", Calo_EventID, nCaloHead/4.);
	logger->info("CsITK: {} expected {}", CsI_EventID, nCsIHead);
	logger->info("HK: {} expected {}", TPointID, nHKHead/4.);
    fout->cd();
    caloTree->Write();
	csiTree->Write();
	hkTree->Write();
    fout->Close();
}

void ComReader::getOutputName(const std::string& prefix,const std::string& filename){
    oname = filename.substr(filename.find_last_of("/")+1);
    oname = oname.substr(0,oname.find_last_of("."));
    oname = prefix+oname+".root";
}

ComReader::~ComReader(){
    file->close();
    // delete file;
    
}
