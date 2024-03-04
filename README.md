# dl-ee-ext

## sel.root
```bash
T_sel->Draw("(E_nu_rec-E_nu_tru)/E_nu_tru>>h(100,-1,1)","")
T_sel->Draw("(E_nu_rec-E_nu_tru)/E_nu_tru>>h(100,-1,1)","match_isFC==1&&ke_info==1")

T_sel->Draw("(E_l_rec-E_l_tru)/E_l_tru>>h(100,-1,1)","","")
T_sel->Draw("(E_l_rec-E_l_tru)/E_l_tru>>h(100,-1,1)","match_isFC==1&&ke_info==1")

T_sel->Draw("(E_nu_rec-E_nu_tru)/E_nu_tru:(E_l_rec-E_l_tru)/E_l_tru>>h(100,-1,1,100,-1,1)","ke_info==0","colz")
```