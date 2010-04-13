// $Id: DateParser.java,v 1.4 2009/11/03 19:16:34 cmzmasek Exp $
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
// WWW: www.phylosoft.org/forester

package org.forester.io.parsers.phyloxml.phylogenydata;

import java.math.BigDecimal;

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.Date;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.util.ForesterUtil;

public class DateParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new DateParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private DateParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhylogenyParserException {
        String unit = "";
        if ( element.isHasAttribute( PhyloXmlMapping.CLADE_DATE_UNIT ) ) {
            unit = element.getAttribute( PhyloXmlMapping.CLADE_DATE_UNIT );
        }
        String val = null;
        String min = null;
        String max = null;
        String desc = "";
        for( int j = 0; j < element.getNumberOfChildElements(); ++j ) {
            final XmlElement e = element.getChildElement( j );
            if ( e.getQualifiedName().equals( PhyloXmlMapping.CLADE_DATE_VALUE ) ) {
                val = e.getValueAsString();
            }
            else if ( e.getQualifiedName().equals( PhyloXmlMapping.CLADE_DATE_MIN ) ) {
                min = e.getValueAsString();
            }
            else if ( e.getQualifiedName().equals( PhyloXmlMapping.CLADE_DATE_MAX ) ) {
                max = e.getValueAsString();
            }
            else if ( e.getQualifiedName().equals( PhyloXmlMapping.CLADE_DATE_DESC ) ) {
                desc = e.getValueAsString();
            }
        }
        BigDecimal val_bd = null;
        BigDecimal min_bd = null;
        BigDecimal max_bd = null;
        if ( !ForesterUtil.isEmpty( val ) ) {
            val_bd = new BigDecimal( val );
        }
        if ( !ForesterUtil.isEmpty( min ) ) {
            min_bd = new BigDecimal( min );
        }
        if ( !ForesterUtil.isEmpty( max ) ) {
            max_bd = new BigDecimal( max );
        }
        return new Date( desc, val_bd, min_bd, max_bd, unit );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
