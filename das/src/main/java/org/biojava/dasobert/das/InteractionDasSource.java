/*
 * 
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA, or see the FSF site: http://www.fsf.org.
 * 
 */
package org.biojava.dasobert.das ;

import org.biojava.dasobert.dasregistry.Das1Source;
import org.biojava.dasobert.dasregistry.DasCoordinateSystem;
import org.biojava.dasobert.dasregistry.DasSource;
import org.biojava.utils.xml.XMLWriter; 
import java.io.IOException;
import java.lang.reflect.Method;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.List;


/**
 * Extends the normal DAS source with some interaction specific things
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 *
 */
public class InteractionDasSource extends Das1Source{
    private boolean isActive ; // if the user included this source or not ..
    private boolean isCompatible;
    private boolean isConverted; // if the query has been converted from e.g. uniprot to entrez
    private boolean registered ; // a flag to trace if source comes from registry or from user config
    public static String DEFAULT_NICKNAME = "ownSource";
    public static Capabilities DEFAULT_CAPABILITY = Capabilities.INTERACTION;
        
    
    /**
     * Basic initilization
     *
     */
    public InteractionDasSource() {
    	super();
        isCompatible = true;
        isConverted = false;
        isActive = true;  // default source is actived and used .
        registered = true; // default true = source comes from registry
        setNickname(DEFAULT_NICKNAME);
        ////String[] caps = new String[1];
        //caps[0] = DEFAULT_CAPABILITY;
       ArrayList<Capabilities> caps=new ArrayList();
       caps.add(DEFAULT_CAPABILITY);
        setCapabilities(caps);
    }
    
    public static InteractionDasSource fromDas1Source(Das1Source orig){
    	
    	InteractionDasSource ds = new InteractionDasSource();
    	
    	try {
    	 Class c = Class.forName("org.biojava.dasobert.dasregistry.Das1Source");
    	 Class ic = Class.forName(" org.biojava.dasobert.das.InteractionDasSource");
         Method[] methods  = c.getMethods();
         
         for (int i = 0; i < methods.length; i++) {
             Method m = methods[i];     
             
             String name = m.getName();
             if ( name.substring(0,3).equals("get")) {
                
                 Object o  = m.invoke(orig, new Object[]{});
                                  
                 
                 Method setter = ic.getMethod("set"+name.substring(3,name.length()),new Class[]{} );
                 System.out.println(setter);
                 setter.invoke(ds, o);
             }
             
         }
    	} catch (Exception e){
    		e.printStackTrace();
    	}
         return ds;
    }

   
    /** 
     * a flag if this das source is active
     * @param flag the active flag
     */
    public void setIsActive(boolean flag){ 
    	isActive = flag ; 
    }
    
    /**
     * Specifies whether a source is currentl active i.e. included into the results
     * @return the active flag
     */
    public boolean getIsActive(){ 
    	return isActive; 
    }
    
    /**
     * Checks whether the source is compatible to the coordinate system passed over
     * @return True if the coordinate system matches, false otherwise
     */
    public boolean getIsCompatible(){
    	return this.isCompatible;
    }
    
    public void setIsCompatible(boolean val){
    	this.isCompatible = val;
    }
    
    public boolean getIsConverted(){
    	return this.isConverted;
    }
    
    public void setIsConverted(boolean conv){
    	this.isConverted = conv;
    }
    
    /**
     * Checks whether the source is compatible to the coordinate system passed over
     * @param queryCoordSys The coordinate system of the query
     * @return True if the coordinate system matches, false otherwise
     */
    public boolean getIsCompatible(DasCoordinateSystem queryCoordSys){
    	if (queryCoordSys != null){
    		List<DasCoordinateSystem> cse = this.getCoordinateSystem();
    		for (int i = 0; i < cse.size(); i++){
    			if (cse.get(i).equals(queryCoordSys)){
    				this.isCompatible = true;
    				return true;
    			}
    		}
    	}
    	this.isCompatible = false;
    	return false;
    }

    
    /**
     * 
     * @param flag the registered flag to set
     */
    public void setRegistered(boolean flag){ 
    	registered = flag ; 
    }
    
    /**
     * 
     * @return the registered flag
     */
    public boolean getRegistered(){ 
    	return registered ; 
    }
    
    
    /** 
     * converts a das source to an interaction das source 
     * 
     * @param ds a DasSource to be converted
     * @return a new InteractionDasSource object
     * */
    public static InteractionDasSource  fromDasSource(DasSource ds) {
        InteractionDasSource ids = new InteractionDasSource();
        ids.setUrl(ds.getUrl());
        ids.setAdminemail(ds.getAdminemail());
        ids.setDescription(ds.getDescription());
        ids.setCoordinateSystem(ds.getCoordinateSystem());
        ids.setCapabilities(ds.getCapabilities());
        ids.setRegisterDate(ds.getRegisterDate());
        ids.setLeaseDate(ds.getLeaseDate());
        ids.setNickname(ds.getNickname());
        ids.setId(ds.getId());
        ids.setLabels(ds.getLabels());
        ids.setHelperurl(ds.getHelperurl());
        return ids;
    }
    
    
    /**
     * Simple to string
     * @return the string
     */
    public String toString() {
        String txt = getId()  + " " + getNickname() + " " + getUrl() ;
        return txt;
    }
    
    /** convert to XML.
     * 
     * @param xw an XMLWriter to write to
     * @return XMLWriter returns it again
     * @throws IOException
     *  */
    public XMLWriter toXML(XMLWriter xw)
    throws IOException
    {
        //System.out.println("writing XML of" + getUrl());
        xw.openTag("InteractionDasSource");
        xw.attribute("url",getUrl());
        xw.attribute("nickname",getNickname());
        xw.attribute("adminemail",getAdminemail());
        //xw.attribute("description",getDescription());
        xw.attribute("status",""+isActive);
        xw.attribute("registered",""+registered);
        DateFormat df = new SimpleDateFormat("dd.MM.yyyy"); 
        //DateFormat df = DateFormat.getDateInstance();
        String rds = df.format(getRegisterDate());
        String lds = df.format(getLeaseDate());
        xw.attribute("registerDate",rds);
        xw.attribute("leaseDate",lds);
        
        // description
        xw.openTag("description");
        xw.print(getDescription());
        xw.closeTag("description");
        
        // coordinateSystems
        xw.openTag("coordinateSystems");
        List<DasCoordinateSystem> coords = getCoordinateSystem();
        for (int i = 0; i < coords.size(); i++){
            xw.openTag("coordinateSystem");
            xw.attribute("name",coords.get(i).toString());
            xw.attribute("uniqId",coords.get(i).getUniqueId());
            xw.closeTag("coordinateSystem");
        } 
        xw.closeTag("coordinateSystems");
        
        // capabilities
        xw.openTag("capabilities");
        List<Capabilities> caps = getCapabilities();
        for (int i = 0; i < caps.size(); i++){
            xw.openTag("capability");
            xw.attribute("service",caps.get(i).toString());
            xw.closeTag("capability");
        } 
        xw.closeTag("capabilities");
        
        xw.closeTag("InteractionDasSource");
        return xw ;
    }
    
    
}
