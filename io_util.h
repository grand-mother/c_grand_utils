int write_GRANDevent(GRANDEvent *event);
void grand_root_read_runhdr(TFile *g);
void grand_root_read_runsimdata(TFile *g);
int grand_root_read_event(TFile *g,int ievt, GRANDEvent *event);
