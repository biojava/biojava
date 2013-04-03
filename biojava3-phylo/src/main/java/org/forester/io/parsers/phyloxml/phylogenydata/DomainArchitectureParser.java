// $Id: DomainArchitectureParser.java,v 1.2 2009/11/03 19:16:34 cmzmasek Exp $
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

import org.forester.io.parsers.PhylogenyParserException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.XmlElement;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.ProteinDomain;

public class DomainArchitectureParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new DomainArchitectureParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private DomainArchitectureParser() {
    }

    @Override
    public DomainArchitecture parse( final XmlElement element ) throws PhylogenyParserException {
        final DomainArchitecture architecure = new DomainArchitecture();
        if ( !element.isHasAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_LENGTH ) ) {
            throw new PhylogenyParserException( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_LENGTH
                    + " attribute is required for domain architecture" );
        }
        final String lenght_str = element.getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_LENGTH );
        try {
            architecure.setTotalLength( Integer.parseInt( lenght_str ) );
        }
        catch ( final NumberFormatException e ) {
            throw new PhylogenyParserException( "could not extract domain architecture length from [" + lenght_str
                    + "]: " + e.getMessage() );
        }
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_DOMAIN ) ) {
                architecure.addDomain( ( ProteinDomain ) ProteinDomainParser.getInstance().parse( child_element ) );
            }
        }
        return architecure;
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
