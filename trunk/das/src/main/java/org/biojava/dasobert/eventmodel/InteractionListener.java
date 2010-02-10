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

import org.biojava.dasobert.das.InteractionParameters;

/**
 * Interface describing the potential outcomes of an interaction thread
 * @author Hagen Blankenburg, Max Planck Institute for Informatics
 *
 */
public interface InteractionListener 
extends ObjectListener{
    
	/**
	 * Called if new interactions where found
	 * @param event Interactions and source parameters
	 */
	public void newInteractions(InteractionEvent event);
    
	/**
	 * Called if no interactions were found
	 * @param params Source paramters
	 */
    public void noObjectFound(InteractionParameters params);
   
    /**
     * Called if the results have to be prepared first and the client shoudl return later
     */
    public void comeBackLater();
    
    /**
     * Dunno
     */
    void newObjectRequested(String accessionCode);
}
