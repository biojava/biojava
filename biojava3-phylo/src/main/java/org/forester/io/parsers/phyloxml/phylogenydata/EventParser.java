// $Id: EventParser.java,v 1.3 2009/11/03 19:16:34 cmzmasek Exp $
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
import org.forester.phylogeny.data.Confidence;
import org.forester.phylogeny.data.Event;
import org.forester.phylogeny.data.PhylogenyData;
import org.forester.util.ForesterUtil;

public class EventParser implements PhylogenyDataPhyloXmlParser {

    private static final PhylogenyDataPhyloXmlParser _instance;
    static {
        try {
            _instance = new EventParser();
        }
        catch ( final Throwable e ) {
            throw new RuntimeException( e.getMessage() );
        }
    }

    private EventParser() {
    }

    @Override
    public PhylogenyData parse( final XmlElement element ) throws PhylogenyParserException {
        String type = "";
        Confidence conf = null;
        int duplications = Event.DEFAULT_VALUE;
        int speciations = Event.DEFAULT_VALUE;
        int losses = Event.DEFAULT_VALUE;
        for( int i = 0; i < element.getNumberOfChildElements(); ++i ) {
            final XmlElement child_element = element.getChildElement( i );
            if ( child_element.getQualifiedName().equals( PhyloXmlMapping.EVENT_TYPE ) ) {
                type = child_element.getValueAsString();
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.CONFIDENCE ) ) {
                conf = ( ( Confidence ) ConfidenceParser.getInstance().parse( child_element ) );
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.EVENT_DUPLICATIONS ) ) {
                duplications = child_element.getValueAsInt();
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.EVENT_SPECIATIONS ) ) {
                speciations = child_element.getValueAsInt();
            }
            else if ( child_element.getQualifiedName().equals( PhyloXmlMapping.EVENT_LOSSES ) ) {
                losses = child_element.getValueAsInt();
            }
        }
        Event event = null;
        if ( ForesterUtil.isEmpty( type ) ) {
            event = new Event( duplications, speciations, losses );
        }
        else {
            try {
                event = new Event( duplications, speciations, losses, type );
            }
            catch ( final Exception e ) {
                throw new PhylogenyParserException( "problem with " + element.toString() + ": " + e.getMessage() );
            }
        }
        if ( conf != null ) {
            event.setConfidence( conf );
        }
        return event;
    }

    public static PhylogenyDataPhyloXmlParser getInstance() {
        return _instance;
    }
}
