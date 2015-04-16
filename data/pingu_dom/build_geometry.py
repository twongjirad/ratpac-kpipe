import os,sys
from math import sin,cos,pi
import StringIO

# This script generated GDML geometry files for KPIPE. It places the SiPMs parametrically.
# the generated portion is sandwiched between two parts (part1 and part2) below.

part1="""<?xml version="1.0" encoding="UTF-8"?>
<gdml xmlns:gdml="http://cern.ch/2001/Schemas/GDML"	
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd" >

<define>
  <position name="sipm_offset" unit="mm" x="0" y="1.0" z="0"/>
  <rotation name="identity"/>
</define>

<materials>
  <element name="bromine" formula="Br" Z="35"> <atom value="79.904"/> </element>
  <element name="hydrogen" formula="H" Z="1">  <atom value="1.0079"/> </element>
  <element name="nitrogen" formula="N" Z="7">  <atom value="14.0067"/> </element>
  <element name="oxygen" formula="O" Z="8">  <atom value="15.999"/> </element>
  <element name="aluminum" formula="Al" Z="13"> <atom value="26.9815"/>  </element>
  <element name="silicon" formula="Si" Z="14"> <atom value="28.0855"/>  </element>
  <element name="carbon" formula="C" Z="6">  <atom value="12.0107"/>  </element>
  <element name="potassium" formula="K" Z="19"> <atom value="39.0983"/>  </element>
  <element name="chromium" formula="Cr" Z="24"> <atom value="51.9961"/>  </element>
  <element name="iron" formula="Fe" Z="26"> <atom value="55.8450"/>  </element>
  <element name="nickel" formula="Ni" Z="28"> <atom value="58.6934"/>  </element>
  <element name="calcium" formula="Ca" Z="20"> <atom value="40.078"/>   </element>
  <element name="magnesium" formula="Mg" Z="12"> <atom value="24.305"/>   </element>
  <element name="sodium" formula="Na" Z="11"> <atom value="22.99"/>    </element>
  <element name="titanium" formula="Ti" Z="22"> <atom value="47.867"/>   </element>
  <element name="argon" formula="Ar" Z="18"> <atom value="39.9480"/>  </element>
  
  <material Z="1" formula=" " name="cryostat_vacuum">
    <D value="1.e-25" unit="g/cm3"/>
    <atom value="1.0079"/>
  </material>

  <material name="stainless_steel" formula="stainless_steel">
    <D value="7.9300" unit="g/cm3"/>
    <fraction n="0.0010" ref="carbon"/>
    <fraction n="0.1792" ref="chromium"/>
    <fraction n="0.7298" ref="iron"/>
    <fraction n="0.0900" ref="nickel"/>
  </material>
  
  <material formula=" " name="air">
    <D value="0.001205" unit="g/cc"/>
    <fraction n="0.781154" ref="nitrogen"/>
    <fraction n="0.209476" ref="oxygen"/>
    <fraction n="0.00937" ref="argon"/>
  </material>
  
  <material formula=" " name="Dirt">
    <D value="1.7" unit="g/cc"/>
    <fraction n="0.438" ref="oxygen"/>
    <fraction n="0.257" ref="silicon"/>
    <fraction n="0.222" ref="sodium"/>
    <fraction n="0.049" ref="aluminum"/>
    <fraction n="0.019" ref="iron"/>
    <fraction n="0.015" ref="potassium"/>
  </material>
  
  <material formula=" " name="mineral_oil">
    <D value="0.77" unit="g/cc"/>
    <fraction n="0.8563" ref="carbon"/>
    <fraction n="0.1437" ref="hydrogen"/>
  </material>

  <material formula=" " name="pseudocumene">
    <D value="0.8758" unit="g/cc"/>
    <fraction n="0.8994" ref="carbon"/>
    <fraction n="0.1006" ref="hydrogen"/>
  </material>
  
  <material formula=" " name="ppo">
    <D value="1.06" unit="g/cc"/>
    <fraction n="0.8142" ref="carbon"/>
    <fraction n="0.0501" ref="hydrogen"/>
    <fraction n="0.0633" ref="nitrogen"/>
    <fraction n="0.0723" ref="oxygen"/>
  </material>
  
  <material formula=" " name="scintillator">
    <D value="0.78" unit="g/cc"/>
    <fraction n="0.996984" ref="mineral_oil"/>
    <fraction n="0.001919" ref="pseudocumene"/>
    <fraction n="0.001097" ref="ppo"/>
  </material>

  <material formula=" " name="chip_silicon">
    <D value="2.3" unit="g/cc"/>
    <fraction n="1.0" ref="silicon"/>
  </material>
</materials>

<solids>
  <box name="world"
       lunit="cm"
       x="100.0"
       y="100.0"
       z="100.0" />
  <sphere name="s1" lunit="cm" aunit="deg"
	  rmax="10.0" deltaphi="360" deltatheta="180"/>

  <!-- SIPM SOLIDS -->
  <box name="sipm_active"
       lunit="mm"
       x="6"
       y="2"
       z="6" />
  <box name="sipm_active_sub"
       lunit="mm"
       x="6"
       y="2"
       z="6" />
  <box name="sipm_package"
       lunit="mm"
       x="10"
       y="4"
       z="10" />
  <subtraction name="sipm_inactive">
       <first ref="sipm_package"/>
       <second ref="sipm_active_sub"/>
       <positionref ref="sipm_offset"/>
       <rotationref ref="identity"/>
  </subtraction>

</solids>

<structure>
  <!-- building the world inside-out -->

  <!-- define SiPM volumes -->
  <volume name="volActiveSiPM">
    <materialref ref="scintillator"/>
    <solidref ref="sipm_active"/>
  </volume>
  <volume name="volInactiveSiPM">
    <materialref ref="chip_silicon"/>
    <solidref ref="sipm_inactive"/>
  </volume>

  <!-- below is procedurally generated -->
"""

part2 = """
  <!-- above s procedurally generated -->
    <volume name="volWorld">
    <solidref ref="world"/>
    <materialref ref="cryostat_vacuum"/>
    <physvol name="pvS1">
      <volumeref ref="volS1"/>
      <position name="posS1" x="0" y="0" z="0"/>
    </physvol>
  </volume>
</structure>

<setup name="Default" version="1.0">
  <world ref="volWorld" />
</setup>

</gdml>
"""

def get_sipm_pos( nsipms, sipm_radius, pmtinfofile ):
    """generates PMTINFO.db and rturns dict with sipm location and rotations
    
    This is needed because we hacked RAT.  It stil thinks we are using PMTs,
    so we have to give it pseudo-data. Fixing this is on a to-do list somewhere.
    We also take the opportunity to build a look up table of channel and position.

    """
    pmtinfo = StringIO.StringIO()
    print >> pmtinfo, "{"
    print >> pmtinfo, "  name:\"PMTINFO\","
    print >> pmtinfo, "  valid_begin: [0,0],"
    print >> pmtinfo, "  valid_end: [0,0],"
    xposlist = "["
    yposlist = "["
    zposlist = "["
    typelist = "["

    sipm_dict = {}
    if nsipms==1:
        sipm_info = { "x1":0.0, "y1":0.0, "z1":-(sipm_radius-0.1),
                      "x2":0.0, "y2":0.0, "z2":-sipm_radius,
                      "rotx":-90*pi/180.0, "roty":0, "rotz":0 }
        sipm_dict[0] = sipm_info


    for isipm in xrange(0,nsipms):
        xposlist += "%.2f,"%(sipm_dict[isipm]['x1'])
        yposlist += "%.2f,"%(sipm_dict[isipm]['y1'])
        zposlist += "%.2f,"%(sipm_dict[isipm]['z1'])
        typelist += "1,"
        
        xposlist += "]"
    yposlist += "]"
    zposlist += "]"
    typelist += "]"
    print >> pmtinfo, " x:",xposlist,","
    print >> pmtinfo, " y:",yposlist,","
    print >> pmtinfo, " z:",zposlist,","
    print >> pmtinfo," type:",typelist,","
    print >> pmtinfo,"}"
    
    f = open( pmtinfofile,'w')
    print >> f, pmtinfo.getvalue()
        
    return sipm_dict
    

def generate_gdml_file( gdml_filename, pmtinfo_filename, nsipms, sipm_radius=8.0 ):
    """ Generates pingu_dom.gdml file. Populates detector with spatially arranged sipms

    For each SiPM, we lay down two components, the active and inactive SiPMs.
    We do this for Chroma. All units in cm.
    """
    ip_sipmdict = get_sipm_pos( nsipms, sipm_radius, pmtinfo_filename )
    targetvol = "  <volume name=\"volS1\">\n"
    targetvol+= "    <materialref ref=\"scintillator\"/>\n"
    targetvol+= "    <solidref ref=\"s1\"/>\n"
    keys = ip_sipmdict.keys()
    keys.sort()
    for key in keys:
        sipminfo = ip_sipmdict[key]
        # active
        targetvol+="    <physvol name=\"SiPM%d\">\n"%(key)
        targetvol+="      <volumeref ref=\"volActiveSiPM\"/>\n"
        targetvol+="      <rotation name=\"rotVolActiveSiPM%d\" x=\"%.2f\" y=\"%.2f\" z=\"%.2f\"/>\n"%(key,  
                                                                                                       sipminfo['rotx'],
                                                                                                       sipminfo['roty'],
                                                                                                       sipminfo['rotz'] )
        targetvol+="      <position name=\"posVolActiveSiPM%d\" unit=\"cm\" x=\"%.4f\" y=\"%.4f\" z=\"%.4f\"/>\n"%(key, 
                                                                                                                   sipminfo['x1'],
                                                                                                                   sipminfo['y1'],
                                                                                                                   sipminfo['z1'] )
        targetvol+="    </physvol>\n"
        # inactive
        targetvol+="    <physvol name=\"InactiveSiPM%d\">\n"%(key)
        targetvol+="      <volumeref ref=\"volInactiveSiPM\"/>\n"
        targetvol+="      <rotation name=\"rotVolInactiveSiPM%d\" x=\"%.2f\" y=\"%.2f\" z=\"%.2f\"/>\n"%(key, 
                                                                                                         sipminfo['rotx'],
                                                                                                         sipminfo['roty'],
                                                                                                         sipminfo['rotz'] )
        targetvol+="      <position name=\"posVolInactiveSiPM%d\" unit=\"cm\" x=\"%.4f\" y=\"%.4f\" z=\"%.4f\"/>\n"%(key, 
                                                                                                                     sipminfo['x2'],
                                                                                                                     sipminfo['y2'],
                                                                                                                     sipminfo['z2'] )
        targetvol+="    </physvol>\n"
        targetvol+="  </volume>\n"
            
    fgdml = open( gdml_filename, 'w' )
    print >> fgdml, part1+"\n"+targetvol+"\n"+part2+"\n"
    fgdml.close()

    return part1+"\n"+targetvol+"\n"+part2+"\n"


if __name__=="__main__":
    generate_gdml_file( "pingu_dom_1sipm.gdml", "PMTINFO.ratdb", 1 )
        

    
