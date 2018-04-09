import os
#from amino_acids_sorted_native_frequency import plotting_amino_acid_frequncy
#------------------------------------
#Store in secondary structures in a counter list
#------------------------------------
secondary_structure=[]
secondary_structure_list=['Turn,Coil',
'Coil,Bridge',
'Coil,Helix',
'Strand,Coil',
'Helix,Bridge',
'Bridge,Turn',
'Strand,Bridge',	
'Turn,Turn',
'Turn,Helix',
'Strand,Helix',
'Bridge,Bridge',
'Strand,Strand',
'Strand,Turn',
'Helix,Helix',
'Coil,Coil']


#Make one large dictionary,secondary structure combination as key:
cys1_secondary_structure_dict_phi={}
cys1_secondary_structure_dict_psi={}

cys2_secondary_structure_dict_phi={}
cys2_secondary_structure_dict_psi={}

get=open('./../'+'angles_secondary_structure.txt','r')
for line in get:
	line=line.strip('\n')
	lines=line.split(',')
	#---
	#Have to store in lists based on secondary structure
	#Will assign the line based to a dictionary:
	
	#Secondary structure list will determine the order of dictionaries
	cys1_ss=lines[-2:]
	cys2_ss=lines[-1:]+ lines[-2:-1]
	cys1_ss=",".join(cys1_ss)
	cys2_ss=",".join(cys2_ss)



	#Phi and Psi angles
	phi    =float(lines[22])
	psi    =float(lines[23])
	phi_b  =float(lines[34])
	psi_b  =float(lines[35])
	phi_a  =float(lines[44])
	psi_a  =float(lines[45])
 
	phi_x  =float(lines[25])
	psi_x  =float(lines[26])
	phi_x_b=float(lines[54])
	psi_x_b=float(lines[55])
	phi_x_a=float(lines[64])
	psi_x_a=float(lines[65])


	if cys1_ss in secondary_structure_list:
		#printprint 'YES',cys1_ss
		try:
			cys1_secondary_structure_dict_phi[cys1_ss].append(phi)
			cys1_secondary_structure_dict_phi[cys1_ss].append(phi)

			#Append Cys2 to 'Cys2 dict'
			cys2_secondary_structure_dict_phi[cys1_ss].append(phi_x)
			cys2_secondary_structure_dict_psi[cys1_ss].append(psi_x)
		except KeyError:
			cys1_secondary_structure_dict_phi[cys1_ss]=[phi]
			cys1_secondary_structure_dict_psi[cys1_ss]=[psi]

			cys2_secondary_structure_dict_phi[cys1_ss]=[phi_x]
			cys2_secondary_structure_dict_psi[cys1_ss]=[psi_x]

	
	#If CYS1_ss not in secondary strcutre list then the order of secondary structure is reveresed. i.e (Strand Helix instetad of Helix Strand) and then Cys1 and Cys2 are also swapped
	if cys1_ss not in secondary_structure_list:
		try:
			cys1_secondary_structure_dict_phi[cys2_ss].append(phi_x)
			cys1_secondary_structure_dict_phi[cys2_ss].append(phi_x)

			#Append Cys2 to 'Cys2 dict'
			cys2_secondary_structure_dict_phi[cys2_ss].append(phi)
			cys2_secondary_structure_dict_psi[cys2_ss].append(psi)
		except KeyError:
			cys1_secondary_structure_dict_phi[cys2_ss]=[phi_x]
			cys1_secondary_structure_dict_psi[cys2_ss]=[psi_x]

			cys2_secondary_structure_dict_phi[cys2_ss]=[phi]
			cys2_secondary_structure_dict_psi[cys2_ss]=[psi]


for key in cys1_secondary_structure_dict_phi:
	print key
	#plotting_amino_acid_frequncy(secondary_structure_dict[key],key)


xdr=range(-180,181)
ydr=range(-180,181)
def make_plot(xname, yname, xax=False, yax=False):
    mbl = 40 if yax else 0
    mbb = 40 if xax else 0


    plot = Plot(
        x_range=xdr, y_range=ydr, background_fill_color="#efe8e2",
        border_fill_color='white', plot_width=200 + mbl, plot_height=200 + mbb,
        min_border_left=2+mbl, min_border_right=2, min_border_top=2, min_border_bottom=2+mbb)

    circle = Circle(x=xname, y=yname, fill_color="color", fill_alpha=0.2, size=4, line_color="color")
    r = plot.add_glyph(source, circle)

    xdr.renderers.append(r)
    ydr.renderers.append(r)

    xticker = BasicTicker()
    if xax:
        xaxis = LinearAxis()
        plot.add_layout(xaxis, 'below')
        xticker = xaxis.ticker
    plot.add_layout(Grid(dimension=0, ticker=xticker))

    yticker = BasicTicker()
    if yax:
        yaxis = LinearAxis()
        plot.add_layout(yaxis, 'left')
        yticker = yaxis.ticker
    plot.add_layout(Grid(dimension=1, ticker=yticker))

    plot.add_tools(PanTool(), WheelZoomTool())

    return plot