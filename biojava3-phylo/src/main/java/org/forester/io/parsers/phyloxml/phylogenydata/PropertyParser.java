// $Id: PropertyParser.java,v 1.3 2009/11/03 19:16:34 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
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
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/

package org.forester.io.parsers.phyloxml.phylogenydata;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.phylogeny.data.Property;
import org.forester.phylogeny.data.Property.AppliesTo;
import org.forester.util.ForesterUtil;

public class PropertyParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new PropertyParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private PropertyParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhylogenyParserException {
        String ref = "";
        String value = "";
        String unit = "";
        String datatype = "";
        String applies_to_str = "";
        String id_ref = "";
        if ( element.isHasAttribute( PhyloXmlMapping.PROPERTY_REF ) ) {
            ref = element.getAttribute( PhyloXmlMapping.PROPERTY_REF );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.PROPERTY_UNIT ) ) {
            unit = element.getAttribute( PhyloXmlMapping.PROPERTY_UNIT );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.PROPERTY_DATATYPE ) ) {
            datatype = element.getAttribute( PhyloXmlMapping.PROPERTY_DATATYPE );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.PROPERTY_APPLIES_TO ) ) {
            applies_to_str = element.getAttribute( PhyloXmlMapping.PROPERTY_APPLIES_TO );
        }
        if ( element.isHasAttribute( PhyloXmlMapping.ID_REF ) ) {
            id_ref = element.getAttribute( PhyloXmlMapping.ID_REF );
        }
        if ( !ForesterUtil.isEmpty( element.getValueAsString() ) ) {
            value = element.getValueAsString();
        }
        AppliesTo applies_to = AppliesTo.OTHER;
        if ( applies_to_str.equals( AppliesTo.NODE.toString() ) ) {
            applies_to = AppliesTo.NODE;
        }
        else if ( applies_to_str.equals( AppliesTo.PARENT_BRANCH.toString() ) ) {
            applies_to = AppliesTo.PARENT_BRANCH;
        }
        else if ( applies_to_str.equals( AppliesTo.CLADE.toString() ) ) {
            applies_to = AppliesTo.CLADE;
        }
        else if ( applies_to_str.equals( AppliesTo.ANNOTATION.toString() ) ) {
            applies_to = AppliesTo.ANNOTATION;
        }
        else if ( applies_to_str.equals( AppliesTo.PHYLOGENY.toString() ) ) {
            applies_to = AppliesTo.PHYLOGENY;
        }
        return new Property( ref, value, unit, datatype, applies_to, id_ref );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
