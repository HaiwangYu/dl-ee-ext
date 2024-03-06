# dl-ee-ext

## sel.root
```bash
T_sel->Draw("(E_nu_rec-E_nu_tru)/E_nu_tru>>h(100,-1,1)","")
T_sel->Draw("(E_nu_rec-E_nu_tru)/E_nu_tru>>h(100,-1,1)","match_isFC==1&&ke_info==1")

T_sel->Draw("(E_l_rec-E_l_tru)/E_l_tru>>h(100,-1,1)","","")
T_sel->Draw("(E_l_rec-E_l_tru)/E_l_tru>>h(100,-1,1)","match_isFC==1&&ke_info==1")

T_sel->Draw("(E_nu_rec-E_nu_tru)/E_nu_tru:(E_l_rec-E_l_tru)/E_l_tru>>h(100,-1,1,100,-1,1)","ke_info==0","colz")
```

```bash
T_corr->Draw("(E_nu_rec-E_nu_tru)/E_nu_tru>>h(100,-1,1)","match_isFC==1")
T_corr->Draw("(E_nu_rec_corr-E_nu_tru)/E_nu_tru>>h(100,-1,1)","match_isFC==1")


T_corr->Draw("(E_nu_rec_corr-E_nu_tru)/E_nu_tru:E_nu_tru>>h(50,0,5,100,-1,1)","match_isFC==1","colz")

T_corr->Draw("E_nu_tru-E_nu_rec:E_nu_rec>>h(50,0,5,100,-5,5)","match_isFC==1","colz")
T_corr->Draw("E_nu_tru-E_nu_rec_corr:E_nu_rec_corr>>h(50,0,5,100,-5,5)","match_isFC==1","colz")

T_corr->Draw("(E_nu_rec_corr-E_nu_tru)/E_nu_tru:E_nu_rec_corr>>h(50,0,5,100,-1,1)","match_isFC==1","colz")
```