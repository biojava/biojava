// $Id: BasicPhylogeniesWebserviceClient.java,v 1.3 2008/10/05 09:11:50 cmzmasek
// Exp $
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

import org.forester.archaeopteryx.webservices.WebservicesManager.WsPhylogenyFormat;
import org.forester.util.ForesterUtil.PhylogenyNodeField;

public class BasicPhylogeniesWebserviceClient implements PhylogeniesWebserviceClient {

    private final String             _desc;
    private final String             _instructions;
    private final String             _menu_name;
    private final String             _name;
    private final WsPhylogenyFormat  _format;
    private final String             _url;
    private final boolean            _integer;
    private final PhylogenyNodeField _node_field;
    private final Object             _proc_inst;
    private final String             _ref;

    public BasicPhylogeniesWebserviceClient( final String name,
                                             final String menu_name,
                                             final String desc,
                                             final String instructions,
                                             final WsPhylogenyFormat format,
                                             final PhylogenyNodeField node_field,
                                             final String url,
                                             final boolean integer,
                                             final String ref,
                                             final Object proc_inst ) {
        super();
        _desc = desc;
        _instructions = instructions;
        _menu_name = menu_name;
        _name = name;
        _format = format;
        _node_field = node_field;
        _url = url;
        _integer = integer;
        _ref = ref;
        _proc_inst = proc_inst;
    }

    @Override
    public String getDescription() {
        return _desc;
    }

    @Override
    public String getInstructions() {
        return _instructions;
    }

    @Override
    public String getMenuName() {
        return _menu_name;
    }

    @Override
    public String getName() {
        return _name;
    }

    @Override
    public PhylogenyNodeField getNodeField() {
        return _node_field;
    }

    @Override
    public Object getProcessingInstructions() {
        return _proc_inst;
    }

    @Override
    public String getReference() {
        return _ref;
    }

    @Override
    public WsPhylogenyFormat getReturnFormat() {
        return _format;
    }

    @Override
    public String getUrl() {
        return _url;
    }

    @Override
    public boolean isQueryInteger() {
        return _integer;
    }
}
