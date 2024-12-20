#plot pcolor for MAE skill score
fig = plt.figure()
plt.pcolor(MAE_SKILL_MN, cmap=colormap, vmin=cblim_mae*-1, vmax=cblim_mae)
plt.xticks(np.array(range(len(variables)))+0.5,variables, fontsize = labelsize+10)
plt.yticks(np.array(range(len(xlabel)))+0.5,xlabel, fontsize = labelsize+10)
cb = plt.colorbar()
cb.set_label(label='MAE skill score (%)',weight='bold',fontsize=labelsize+10)
cb.ax.tick_params(labelsize=labelsize+8)
fig.savefig('./boxplots/'+stationtype+'/pcolor_MAE_'+valmode+'_'+aggval+'_'+stationtype+extension, dpi=300)
plt.close(fig)

#plot pcolor for HIMAE_obs skill score
fig = plt.figure()
plt.pcolor(HIMAE_obs_SKILL_MN, cmap=colormap, vmin=cblim_himae_obs*-1, vmax=cblim_himae_obs)
plt.xticks(np.array(range(len(variables)))+0.5,variables, fontsize = labelsize+10)
plt.yticks(np.array(range(len(xlabel)))+0.5,xlabel, fontsize = labelsize+10)
cb = plt.colorbar()
cb.set_label(label='HIMAE_obs skill score (%)',weight='bold',fontsize=labelsize+10)
cb.ax.tick_params(labelsize=labelsize+8)
fig.savefig('./boxplots/'+stationtype+'/pcolor_HIMAE_obs_'+valmode+'_'+aggval+'_'+stationtype+extension, dpi=300)
plt.close(fig)

#plot pcolor for LOMAE_obs skill score
fig = plt.figure()
plt.pcolor(LOMAE_obs_SKILL_MN, cmap=colormap, vmin=cblim_himae_obs*-1, vmax=cblim_himae_obs)
plt.xticks(np.array(range(len(variables)))+0.5,variables, fontsize = labelsize+10)
plt.yticks(np.array(range(len(xlabel)))+0.5,xlabel, fontsize = labelsize+10)
cb = plt.colorbar()
cb.set_label(label='LOMAE_obs skill score (%)',weight='bold',fontsize=labelsize+10)
cb.ax.tick_params(labelsize=labelsize+8)
fig.savefig('./boxplots/'+stationtype+'/pcolor_LOMAE_obs_'+valmode+'_'+aggval+'_'+stationtype+extension, dpi=300)
plt.close(fig)

#plot pcolor for HIMAE_mod skill score
fig = plt.figure()
plt.pcolor(HIMAE_mod_SKILL_MN, cmap=colormap, vmin=cblim_himae_mod*-1, vmax=cblim_himae_mod)
plt.xticks(np.array(range(len(variables)))+0.5,variables, fontsize = labelsize+10)
plt.yticks(np.array(range(len(xlabel)))+0.5,xlabel, fontsize = labelsize+10)
cb = plt.colorbar()
cb.set_label(label='HIMAE_mod skill score (%)',weight='bold',fontsize=labelsize+10)
cb.ax.tick_params(labelsize=labelsize+8)
fig.savefig('./boxplots/'+stationtype+'/pcolor_HIMAE_mod_'+valmode+'_'+aggval+'_'+stationtype+extension, dpi=300)
plt.close(fig)

#plot pcolor for LOMAE_mod skill score
fig = plt.figure()
plt.pcolor(LOMAE_mod_SKILL_MN, cmap=colormap, vmin=cblim_himae_mod*-1, vmax=cblim_himae_mod)
plt.xticks(np.array(range(len(variables)))+0.5,variables, fontsize = labelsize+10)
plt.yticks(np.array(range(len(xlabel)))+0.5,xlabel, fontsize = labelsize+10)
cb = plt.colorbar()
cb.set_label(label='LOMAE_mod skill score (%)',weight='bold',fontsize=labelsize+10)
cb.ax.tick_params(labelsize=labelsize+8)
fig.savefig('./boxplots/'+stationtype+'/pcolor_LOMAE_mod_'+valmode+'_'+aggval+'_'+stationtype+extension, dpi=300)
plt.close(fig)


#calculate performance gain for last variable
horup_vert10 = (np.mean(np.median(MAE[[4,6],:],axis=1)) - np.mean(np.median(MAE[[0,2],:],axis=1)))/ np.mean(np.median(MAE[[0,2],:],axis=1))*100
horup_vert20 = (np.mean(np.median(MAE[[5,7],:],axis=1)) - np.mean(np.median(MAE[[1,3],:],axis=1)))/ np.mean(np.median(MAE[[1,3],:],axis=1))*100

vertup_horc = (np.mean(np.median(MAE[[1,3],:],axis=1)) - np.mean(np.median(MAE[[0,2],:],axis=1)))/ np.mean(np.median(MAE[[0,2],:],axis=1))*100
vertup_horf = (np.mean(np.median(MAE[[5,7],:],axis=1)) - np.mean(np.median(MAE[[4,6],:],axis=1)))/ np.mean(np.median(MAE[[4,6],:],axis=1))*100

tosaprc_vert10 = (np.mean(np.median(MAE[[2,6],:],axis=1)) - np.mean(np.median(MAE[[0,4],:],axis=1)))/ np.mean(np.median(MAE[[0,4],:],axis=1))*100
tosaprc_vert20 = (np.mean(np.median(MAE[[3,7],:],axis=1)) - np.mean(np.median(MAE[[1,5],:],axis=1)))/ np.mean(np.median(MAE[[1,5],:],axis=1))*100

tosaprc_horc = (np.mean(np.median(MAE[[2,3],:],axis=1)) - np.mean(np.median(MAE[[0,1],:],axis=1)))/ np.mean(np.median(MAE[[0,1],:],axis=1))*100
tosaprc_horf = (np.mean(np.median(MAE[[6,7],:],axis=1)) - np.mean(np.median(MAE[[5,6],:],axis=1)))/ np.mean(np.median(MAE[[5,6],:],axis=1))*100

MAEmed = np.median(MAE,axis=1)
perdev = (MAEmed[1:]-MAEmed[0]) / MAEmed[0]*100
print(perdev)
