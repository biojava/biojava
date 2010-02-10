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
package org.biojava.dasobert.eventmodel;

import de.mpg.mpiinf.ag3.dasmi.model.Interaction;
import org.biojava.dasobert.das.InteractionParameters;


/**
 * Interaction event, containing the resutls of an interaction thread
 * @author Hagen Blankenburg
 *
 */
public class InteractionEvent extends AbstractDasEvent {
    Interaction[] interactions = null;
    InteractionParameters params = null;
    
    /**
     * Creates a new InteractionEvent object with the parameters used and the interactions found
     * @param params
     * @param interactions
     */
    public InteractionEvent(InteractionParameters params, Interaction[] interactions){
        super();
        this.interactions = interactions;
        this.params = params;
    }
    
    /**
     * 
     * @return The InteractionParameters used
     */
    public InteractionParameters getParams(){
        return this.params;
    }
    
    /**
     * 
     * @return The interactions returned from a server
     */
    public Interaction[] getInteractions(){
        return interactions;
    }

}
