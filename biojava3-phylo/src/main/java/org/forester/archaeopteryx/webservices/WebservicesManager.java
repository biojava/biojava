// $Id: WebservicesManager.java,v 1.2 2009/10/28 19:11:22 cmzmasek Exp $
// forester -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2010 Christian M. Zmasek
// Copyright (C) 2008-2010 Burnham Institute for Medical Research
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
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org

package org.forester.archaeopteryx.webservices;

import java.util.ArrayList;
import java.util.List;

public final class WebservicesManager {

    private static WebservicesManager               _instance;
    private final List<PhylogeniesWebserviceClient> _clients;

    private WebservicesManager() {
        _clients = new ArrayList<PhylogeniesWebserviceClient>();
        _clients.addAll( WebserviceUtil.createDefaultClients() );
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        throw new CloneNotSupportedException();
    }

    public PhylogeniesWebserviceClient getAvailablePhylogeniesWebserviceClient( final int i ) {
        return getAvailablePhylogeniesWebserviceClients().get( i );
    }

    public List<PhylogeniesWebserviceClient> getAvailablePhylogeniesWebserviceClients() {
        return _clients;
    }

    public static WebservicesManager getInstance() {
        if ( _instance == null ) {
            _instance = new WebservicesManager();
        }
        return _instance;
    }

    public enum WsPhylogenyFormat {
        NH, NHX, NEXUS, TOL_XML_RESPONSE, PHYLOXML, NH_EXTRACT_TAXONOMY, PFAM
    }
}
