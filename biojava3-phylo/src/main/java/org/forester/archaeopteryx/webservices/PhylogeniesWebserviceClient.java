// $Id: PhylogeniesWebserviceClient.java,v 1.3 2008/10/05 09:11:50 cmzmasek Exp
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2010 Christian M. Zmasek
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek at yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx.webservices;

import org.forester.archaeopteryx.webservices.WebservicesManager.WsPhylogenyFormat;
import org.forester.util.ForesterUtil.PhylogenyNodeField;

/*
 * Webservices which return phylogenies.
 */
public interface PhylogeniesWebserviceClient {

    public final static String QUERY_PLACEHOLDER = "__query__";

    /**
     * A short description of the webservice (~20 characters).
     *  
     * @return a short description of the webservice (~20 characters)
     */
    public String getDescription();

    /**
     * Instructions (and examples) on how to use the webservice.
     * 
     * @return instructions (and examples) on how to use the webservice
     */
    public String getInstructions();

    /**
     * A name/description which can appear on a menu.
     * 
     * @return A name/description which can appear on a menu
     */
    public String getMenuName();

    /**
     * The name of the webservice.
     * 
     * 
     * @return the name of the webservice
     */
    public String getName();

    /**
     * The node data field in which to place node names from simple unannotated formats
     * (such as Newick). Null means avoiding any such postprocessing.  
     * 
     * @return the field code
     */
    public PhylogenyNodeField getNodeField();

    /**
     * This is used to indicate any kind of special processing.
     * 
     * 
     * @return a reference
     */
    public Object getProcessingInstructions();

    /**
     * To get a type of reference for the webservice (an URL or citation, for example).
     * 
     * 
     * @return a reference
     */
    public String getReference();

    /**
     * The expected format of the response.
     * 
     * @return the expected format of the response
     */
    public WsPhylogenyFormat getReturnFormat();

    /**
     * Use QUERY_PLACEHOLDER to indicate position of query variable.
     * 
     * @return the URL
     */
    public String getUrl();

    /**
     * Is the query a number?
     * 
     * 
     * @return
     */
    public boolean isQueryInteger();
}
