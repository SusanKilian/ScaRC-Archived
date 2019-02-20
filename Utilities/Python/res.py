#me -*- coding: utf-8 -*-

#####################################################################################
# Initial import
#####################################################################################
import os.path

import sys
import math
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt


def  read_csv(solver, kchid, kite_total, kism, knl, kmethod, ksolver, 
              kite, kite_cg, kite_mg, kite_sm, 
              dres, dres_cg, dres_mg, dres_sm, 
              derr, derr_cg, derr_mg, derr_sm) :


    kchid0 = []
    kite_total0 = []
    
    kism0 = []
    knl0 = []
    kmethod0 = []
    ksolver0 = []
    
    kite0 = []
    kite_cg0 = []
    kite_mg0 = []
    kite_sm0 = []
    
    dres0 = []
    dres_cg0 = []
    dres_mg0 = []
    dres_sm0 = []
    
    derr0 = []
    derr_cg0 = []
    derr_mg0 = []
    derr_sm0 = []

    name_chid= "res/%s_residual.csv" % (solver)
    chid0 = name_chid

    if not os.path.isfile(name_chid): 
       print 'Not found:  %s' %name_chid
       return 0

    print 'Processing: %s' %name_chid
    chid_in  = open(name_chid ,'r')
    chid_input  = chid_in.readlines()
    chid_in.close()

    num=0
    start = 1

    for line in chid_input:

        num=num+1

        line, null = line.split ("\n")
        values = line.split (",")

        if (start > 0):
            start -= 1
            continue

        ite_total  = int(values[0])
        ite        = int(values[1])
        ite_cg     = int(values[2])
        ite_mg     = int(values[3])
        ite_sm = int(values[4])
        ism        = int(values[5])
        nl         = int(values[6])
        method     = int(values[7])
        solver     = int(values[8])
        res        = float(values[9])
        err        = float(values[10])

        #if (l<3): r_min = min(1000.0, r)

        kchid0.append(chid)
        kite_total0.append(ite_total)
        kite0.append(ite)
        if (ism == 6):
           kism0.append(1.25)
        elif (ism == 7):
           kism0.append(1.75)
        else:
           kism0.append(1.5)
        knl0.append(nl)
        kmethod0.append(method)
        ksolver0.append(solver)
        dres0.append(res)
        derr0.append(err)

        if (ism != 0):
            kite_sm0.append(ite_total)
            dres_sm0.append(res)
            derr_sm0.append(err)
        if (method == 1 and solver == 1):
            kite_cg0.append(ite_cg)
            dres_cg0.append(res)
            derr_cg0.append(err)
        if (method == 2 and solver == 1):
            kite_mg0.append(ite_mg)
            dres_mg0.append(res)
            derr_mg0.append(err)
            if len(kite_mg0)>1: 
               print "%6d   %10.4e   %8.3f    %10.4e   %8.3f  " %(ite_mg, res, dres_mg0[-1]/dres_mg0[-2], err, derr_mg0[-1]/derr_mg0[-2])
            else:
               print "\n"
               print "  ITE      Residual     Ratio        Error       Ratio "
               print "----------------------------------------------------------"
               print "%6d   %10.4e               %10.4e  " %(ite_mg, res, err)

    print "\n"

    kchid.append(chid0)
    kite_total.append(kite_total0)
    kism.append(kism0)
    knl.append(knl0)
    kmethod.append(kmethod0)
    ksolver.append(ksolver0)
    kite.append(kite0)
    kite_cg.append(kite_cg0)
    kite_mg.append(kite_mg0)
    kite_sm.append(kite_sm0)
    dres.append(dres0)
    dres_cg.append(dres_cg0)
    dres_mg.append(dres_mg0)
    dres_sm.append(dres_sm0)
    derr.append(derr0)
    derr_cg.append(derr_cg0)
    derr_mg.append(derr_mg0)
    derr_sm.append(derr_sm0)


 

def plot_csv3(solver, ires, chid, ite, quan1, quan2, quan3, quan4, name, mevery, print_log, xlim):

    fig = plt.figure (facecolor='w')
    ax = fig.add_subplot(111)

    nsim = len(chid)

    legsize = fnt.FontProperties(size=11)

    legend = []

    box = ax.get_position()

    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.89])
    ax.set_xlim(0,xlim)
    #ax.set_ylim(0.001,100)
    ax.plot (ite[ires], quan1[ires], markevery = mevery, marker = 'o', linewidth=0.75, linestyle='-', color = "r")
    ax.plot (ite[ires], quan2[ires], markevery = mevery, marker = 'o', linewidth=0.75, linestyle='-', color = "c")
    ax.legend(legend, prop=legsize, loc="lower center", bbox_to_anchor=(0.50, -0.3), fancybox=True, shadow=True, ncol=5)
    if (print_log): ax.set_yscale("log")
    ax.grid(True)
    legend.append('Residual')

    ax2 = ax.twinx()
    ax2.set_xlim(0,xlim)
    #ax.set_ylim(-1.0,5.0)
    ax2.plot (ite[ires], quan3[ires], markevery = mevery, marker = 's', linewidth=0.75, linestyle='-', color = "b")
    ax2.plot (ite[ires], quan4[ires], markevery = mevery, marker = 's', linewidth=0.75, linestyle='-', color = "g")
    legend.append('Level')


    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(11)

    fig.tight_layout()
    output = "pictures/residuals/%s_%d.png" % (name, ires)
    savefig(output)
    plt.close(fig)

def plot_csv2(solver, ires, chid, ite, res, err, name, mevery, print_log, xlim):

    fig = plt.figure (facecolor='w')
    ax = fig.add_subplot(111)

    nsim = len(chid)

    legsize = fnt.FontProperties(size=11)

    legend = []

    box = ax.get_position()

    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.89])
    ax.set_xlim(0,xlim)
    #ax.set_ylim(0.001,100)
    ax.plot (ite[ires], res[ires], markevery = mevery, marker = 'o', linewidth=0.75, linestyle='-', color = "r")
    ax.plot (ite[ires], err[ires], markevery = mevery, marker = 'o', linewidth=0.75, linestyle='-', color = "b")
    ax.legend(legend, prop=legsize, loc="lower center", bbox_to_anchor=(0.50, -0.3), fancybox=True, shadow=True, ncol=5)
    if (print_log): ax.set_yscale("log")
    ax.grid(True)
    legend.append('Residual')
    legend.append('Error')

    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(11)

    fig.tight_layout()
    output = "pictures/residuals/%s_%d.png" % (name, ires)
    savefig(output)
    plt.close(fig)



def plot_csv1(solver, ires, chid, ite, res, name, mevery, print_log, xlim):

    fig = plt.figure (facecolor='w')
    ax = fig.add_subplot(111)

    nsim = len(chid)

    legsize = fnt.FontProperties(size=11)

    index  = []
    legend = []

    #pstart0=b1
    #pend0  =b2

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.89])

    ax.plot (ite[ires], res[ires], markevery = mevery, marker = 'o', linewidth=0.75, linestyle='-', color = "r")
    legend.append(name)

    ax.legend(legend, prop=legsize, loc="lower center", bbox_to_anchor=(0.50, -0.3), fancybox=True, shadow=True, ncol=5)

    ax.grid(True)
    if (print_log): ax.set_yscale("log")
    #ax.set_title('%s meshes'%(top),fontsize=30)

    #ax.set_xlabel('Time [s]',fontsize=16)
    #ax.set_ylabel(r'%s' %title[0][iquan+1].strip('\"'),fontsize=16)

    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(11)

    output = "pictures/residuals/%s_%d.png" % (name, ires)
    savefig(output)
    plt.close(fig)



def plot_csv0(solver, ires, chid, time, time_cg, time_mg, res, err, cg, mg, name, mevery, print_log, xlim):

    fig = plt.figure (facecolor='w')
    ax = fig.add_subplot(111)

    nsim = len(chid)

    legsize = fnt.FontProperties(size=11)

    legend = []

    #pstart0=b1
    #pend0  =b2

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.89])
    #ax.set_xlim(0,xlim)

    ax.plot (time[ires], res[ires], markevery = mevery, marker = 'o', linewidth=0.75, linestyle='-', color = "r")
    legend.append('Res')

    ax.plot (time[ires], err[ires], markevery = mevery, marker = 's', linewidth=0.75, linestyle='-', color = "b")
    legend.append('Err')

    #ax.plot (time_cg[ires], cg[ires], markevery = mevery, marker = 'o', linewidth=0.75, linestyle='-', color = "m")
    #legend.append('Cg')

    ax.plot (time_mg[ires], mg[ires], markevery = mevery, marker = 'o', linewidth=0.75, linestyle='-', color = "g")
    legend.append('MG')

    ax.legend(legend, prop=legsize, loc="lower center", bbox_to_anchor=(0.50, -0.3), fancybox=True, shadow=True, ncol=5)

    ax.grid(True)
    if (print_log): ax.set_yscale("log")
    #ax.set_title('%s meshes'%(top),fontsize=30)
    #ax.set_ylim(0.001,100)

    #ax.set_xlabel('Time [s]',fontsize=16)
    #ax.set_ylabel(r'%s' %title[0][iquan+1].strip('\"'),fontsize=16)

    labels = ax.get_xticklabels() + ax.get_yticklabels()
    for label in labels:
        label.set_size(11)

    output = "pictures/residuals/%s_%d.png" % (name, ires)
    savefig(output)
    plt.close(fig)


xlim = float(sys.argv[1])
nres = 1
solver = 'SCARC_MAIN_MULTIGRID'

chid = []
kite_total = []

kism = []
knl = []
kmethod = []
ksolver = []

kite = []
kite_cg = []
kite_mg = []
kite_sm = []

dres = []
dres_cg = []
dres_mg = []
dres_sm = []

derr = []
derr_cg = []
derr_mg = []
derr_sm = []

for ires in range(nres):
   read_csv(solver, chid, kite_total, kism, knl, kmethod, ksolver, 
            kite, kite_cg, kite_mg, kite_sm, 
            dres, dres_cg, dres_mg, dres_sm, 
            derr, derr_cg, derr_mg, derr_sm) 

print chid
for ires in range(len(chid)):
   plot_csv1(solver, ires, chid, kite_total, knl, 'level', 1, False, xlim)
   plot_csv2(solver, ires, chid, kite_total, dres   , derr   , 'res_err'  , 10000, True, xlim)
   plot_csv2(solver, ires, chid, kite_mg   , dres_mg, derr_mg, 'mg'       ,     1, True, xlim)
   plot_csv2(solver, ires, chid, kite_sm   , dres_sm, derr_sm, 'smooth'   ,     1, True, xlim)
   plot_csv3(solver, ires, chid, kite_total, dres   , derr, knl    , kism, 'res_level',     1, True, xlim)




