// $Id: NodeData.java,v 1.34 2009/11/09 22:42:47 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

import org.forester.util.ForesterUtil;

public class NodeData implements PhylogenyData {

    private Event            _event;
    private Sequence         _sequence;
    private Identifier       _node_identifier;
    private Taxonomy         _taxonomy;
    private Distribution     _distribution;
    private Date             _date;
    private BinaryCharacters _binary_characters;
    private PropertiesMap    _properties;
    private Reference        _reference;

    public NodeData() {
        //Constructor
    }

    public StringBuffer asSimpleText() {
        throw new UnsupportedOperationException();
    }

    public StringBuffer asText() {
        throw new UnsupportedOperationException();
    }

    public PhylogenyData copy() {
        final NodeData new_pnd = new NodeData();
        if ( isHasSequence() ) {
            new_pnd.setSequence( ( Sequence ) getSequence().copy() );
        }
        if ( isHasEvent() ) {
            new_pnd.setEvent( ( Event ) getEvent().copy() );
        }
        if ( isHasNodeIdentifier() ) {
            new_pnd.setNodeIdentifier( ( Identifier ) getNodeIdentifier().copy() );
        }
        if ( isHasTaxonomy() ) {
            new_pnd.setTaxonomy( ( Taxonomy ) getTaxonomy().copy() );
        }
        if ( isHasBinaryCharacters() ) {
            new_pnd.setBinaryCharacters( ( BinaryCharacters ) getBinaryCharacters().copy() );
        }
        if ( isHasReference() ) {
            new_pnd.setReference( ( Reference ) getReference().copy() );
        }
        if ( isHasDistribution() ) {
            new_pnd.setDistribution( ( Distribution ) getDistribution().copy() );
        }
        if ( isHasDate() ) {
            new_pnd.setDate( ( Date ) getDate().copy() );
        }
        if ( isHasProperties() ) {
            new_pnd.setProperties( ( PropertiesMap ) getProperties().copy() );
        }
        return new_pnd;
    }

    public BinaryCharacters getBinaryCharacters() {
        return _binary_characters;
    }

    public Date getDate() {
        return _date;
    }

    public Distribution getDistribution() {
        return _distribution;
    }

    public Event getEvent() {
        return _event;
    }

    public Identifier getNodeIdentifier() {
        return _node_identifier;
    }

    public PropertiesMap getProperties() {
        return _properties;
    }

    public Reference getReference() {
        return _reference;
    }

    public Sequence getSequence() {
        return _sequence;
    }

    public Taxonomy getTaxonomy() {
        return _taxonomy;
    }

    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public boolean isHasBinaryCharacters() {
        return getBinaryCharacters() != null;
    }

    public boolean isHasDate() {
        return ( getDate() != null )
                && ( !ForesterUtil.isEmpty( getDate().getDesc() ) || !ForesterUtil.isNull( getDate().getMax() )
                        || !ForesterUtil.isNull( getDate().getMin() ) || !ForesterUtil.isNull( getDate().getValue() ) || !ForesterUtil
                        .isEmpty( getDate().getUnit() ) );
    }

    public boolean isHasDistribution() {
        return ( getDistribution() != null )
                && ( !ForesterUtil.isEmpty( getDistribution().getDesc() )
                        || !ForesterUtil.isNull( getDistribution().getLatitude() )
                        || !ForesterUtil.isNull( getDistribution().getLongitude() )
                        || !ForesterUtil.isNull( getDistribution().getAltitude() ) || !ForesterUtil
                        .isEmpty( getDistribution().getGeodeticDatum() ) );
    }

    public boolean isHasEvent() {
        return getEvent() != null;
    }

    public boolean isHasNodeIdentifier() {
        return getNodeIdentifier() != null;
    }

    public boolean isHasProperties() {
        return getProperties() != null;
    }

    public boolean isHasReference() {
        return ( getReference() != null )
                && ( !ForesterUtil.isEmpty( getReference().getDoi() ) || !ForesterUtil.isEmpty( getReference()
                        .getValue() ) );
    }

    public boolean isHasSequence() {
        return getSequence() != null;
    }

    public boolean isHasTaxonomy() {
        return getTaxonomy() != null;
    }

    public void setBinaryCharacters( final BinaryCharacters binary_characters ) {
        _binary_characters = binary_characters;
    }

    public void setDate( final Date date ) {
        _date = date;
    }

    public void setDistribution( final Distribution distribution ) {
        _distribution = distribution;
    }

    public void setEvent( final Event event ) {
        _event = event;
    }

    public void setNodeIdentifier( final Identifier node_identifier ) {
        _node_identifier = node_identifier;
    }

    public void setProperties( final PropertiesMap custom_data ) {
        _properties = custom_data;
    }

    public void setReference( final Reference reference ) {
        _reference = reference;
    }

    public void setSequence( final Sequence sequence ) {
        _sequence = sequence;
    }

    public void setTaxonomy( final Taxonomy taxonomy ) {
        _taxonomy = taxonomy;
    }

    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        if ( isHasNodeIdentifier() ) {
            sb.append( getNodeIdentifier().toNHX() );
        }
        if ( isHasTaxonomy() ) {
            sb.append( getTaxonomy().toNHX() );
        }
        if ( isHasSequence() ) {
            sb.append( getSequence().toNHX() );
        }
        if ( isHasEvent() ) {
            sb.append( getEvent().toNHX() );
        }
        if ( isHasProperties() ) {
            sb.append( getProperties().toNHX() );
        }
        return sb;
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( isHasNodeIdentifier() ) {
            writer.write( ForesterUtil.LINE_SEPARATOR );
            writer.write( indentation );
            //  if ( !org.forester.util.ForesterUtil.isEmpty( getNodeIdentifier().getProvider() ) ) {
            //     PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.NODE_IDENTIFIER, getNodeIdentifier()
            //             .getValue(), PhyloXmlMapping.IDENTIFIER_PROVIDER_ATTR, getNodeIdentifier().getProvider() );
            // }
            // else {
            //     PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.NODE_IDENTIFIER, getNodeIdentifier()
            //             .getValue() );
            // }
        }
        if ( isHasTaxonomy() ) {
            getTaxonomy().toPhyloXML( writer, level, indentation );
        }
        if ( isHasSequence() ) {
            getSequence().toPhyloXML( writer, level, indentation );
        }
        if ( isHasEvent() ) {
            getEvent().toPhyloXML( writer, level, indentation );
        }
        if ( isHasBinaryCharacters() ) {
            getBinaryCharacters().toPhyloXML( writer, level, indentation );
        }
        if ( isHasDistribution() ) {
            getDistribution().toPhyloXML( writer, level, indentation );
        }
        if ( isHasDate() ) {
            getDate().toPhyloXML( writer, level, indentation );
        }
        if ( isHasReference() ) {
            getReference().toPhyloXML( writer, level, indentation );
        }
        if ( isHasProperties() ) {
            getProperties().toPhyloXML( writer, level, indentation.substring( 0, indentation.length() - 3 ) );
        }
    }
}
