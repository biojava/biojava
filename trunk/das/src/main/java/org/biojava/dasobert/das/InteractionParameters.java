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
package org.biojava.dasobert.das;

import org.biojava.dasobert.dasregistry.DasCoordinateSystem;

/** 
 * 
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 *
 */
public class InteractionParameters {
    private String[] queries = null;
    private String[] origQueries = null;
    private DasCoordinateSystem queryCoordinateSystem = null;
    private InteractionDasSource dasSource = null;
    private boolean qmSource = false;
        
    /**
     * Empty constructor
     *
     */
    public InteractionParameters() {
        super();
    }


    /**
     * 
     * @return The source associated with an interaction thread 
     */
    public InteractionDasSource getDasSource() {
        return this.dasSource;
    }

    
    /**
     * 
     * @param dasSource The source to set
     */
    public void setDasSource(InteractionDasSource dasSource) {
        this.dasSource = dasSource;
    }

    /**
     * 
     * @return The query
     */
    public String[] getQueries() {
        return this.queries;
    }

    /**
     * 
     * @param queries the queries to set
     */
    public void setQueries(String[] queries) {
        this.queries = queries;
    }
    
    /**
     * 
     * @return The query
     */
    public String[] getOrigQueries() {
        return this.origQueries;
    }

    /**
     * 
     * @param origQueries the queries to set
     */
    public void setOrigQueries(String[] origQueries) {
        this.origQueries = origQueries;
    }


    /**
     * 
     * @return The coordinate system
     */
    public DasCoordinateSystem getQueryCoordinateSystem() {
        return this.queryCoordinateSystem;
    }


    /**
     * 
     * @param queryCoordinateSystem The coordinate system to set
     */
    public void setQueryCoordinateSystem(DasCoordinateSystem queryCoordinateSystem) {
        this.queryCoordinateSystem = queryCoordinateSystem;
    }
    
    
    /**
     * Spcifies if the source specified in the parameters is a qm source only providigng
     * additional confidence score and not interactions 
     * TODO overthink this ...
     * @param isQm the qm flag to set
     */
    public void setQmSource(boolean isQm){
    	this.qmSource = isQm;
    }
    
    /**
     * 
     * @return the qm flag
     */
    public boolean getQmSource(){
    	return this.qmSource;
    }
}
