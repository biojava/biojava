// $Id: ProteinDomainParser.java,v 1.2 2009/11/03 19:16:34 cmzmasek Exp $
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
import org.forester.phylogeny.data.ProteinDomain;

public class ProteinDomainParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new ProteinDomainParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private ProteinDomainParser() {
    }

    @Override
    public ProteinDomain parse( final XmlElement element ) throws PhylogenyParserException {
        String name = "";
        int f = -1;
        int t = -1;
        double conf = ProteinDomain.CONFIDENCE_DEFAULT;
        String id = ProteinDomain.IDENTIFIER_DEFAULT;
        try {
            f = Integer
                    .parseInt( element.getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_FROM ) );
            t = Integer.parseInt( element.getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_TO ) );
            conf = Double.parseDouble( element
                    .getAttribute( PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_PROT_DOMAIN_CONFIDENCE ) );
            if ( element.isHasAttribute( PhyloXmlMapping.IDENTIFIER ) ) {
                id = element.getAttribute( PhyloXmlMapping.IDENTIFIER );
            }
        }
        catch ( final Exception e ) {
            throw new PhylogenyParserException( "failed to parse element [" + element + "]: " + e.getMessage() );
        }
        name = element.getValueAsString();
        if ( ( f == -1 ) || ( t == -1 ) || ( conf == ProteinDomain.CONFIDENCE_DEFAULT ) ) {
            throw new PhylogenyParserException( "from, to, or confidence attribute not set in: " + element );
        }
        return new ProteinDomain( name, f, t, id, conf );
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
