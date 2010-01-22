// $Id: Event.java,v 1.40 2009/10/26 23:29:39 cmzmasek Exp $
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.StringTokenizer;

import org.forester.io.parsers.nhx.NHXFormatException;
import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class Event implements PhylogenyData {

    public final static int     DEFAULT_VALUE = -1;
    private static final String NHX_SEPARATOR = ">";
    private int                 _duplications;
    private int                 _speciations;
    private int                 _gene_losses;
    private EventType           _event_type;
    private Confidence          _confidence;

    public Event() {
        _duplications = DEFAULT_VALUE;
        _speciations = DEFAULT_VALUE;
        _gene_losses = DEFAULT_VALUE;
        _event_type = EventType.unassigned;
    }

    public Event( final EventType type ) {
        _duplications = DEFAULT_VALUE;
        _speciations = DEFAULT_VALUE;
        _gene_losses = DEFAULT_VALUE;
        _event_type = type;
    }

    public Event( final int duplications, final int speciations, final int gene_losses ) {
        _duplications = duplications;
        _speciations = speciations;
        _gene_losses = gene_losses;
        _event_type = EventType.mixed;
    }

    public Event( final int duplications, final int speciations, final int gene_losses, final String type ) {
        _duplications = duplications;
        _speciations = speciations;
        _gene_losses = gene_losses;
        _event_type = EventType.valueOf( type );
    }

    public Event( final String nhx ) throws NHXFormatException {
        if ( ForesterUtil.isEmpty( nhx ) ) {
            _duplications = DEFAULT_VALUE;
            _speciations = DEFAULT_VALUE;
            _gene_losses = DEFAULT_VALUE;
            _event_type = EventType.unassigned;
        }
        else {
            final StringTokenizer st = new StringTokenizer( nhx, NHX_SEPARATOR );
            if ( st.countTokens() != 4 ) {
                throw new NHXFormatException( "malformed NHX format for event [" + nhx + "]" );
            }
            final String duplications = ( String ) st.nextElement();
            final String speciations = ( String ) st.nextElement();
            final String losses = ( String ) st.nextElement();
            final String event_type = ( String ) st.nextElement();
            int d = 0;
            int s = 0;
            int l = 0;
            try {
                d = Integer.parseInt( duplications );
                s = Integer.parseInt( speciations );
                l = Integer.parseInt( losses );
                _duplications = d;
                _speciations = s;
                _gene_losses = l;
                _event_type = EventType.valueOf( event_type );
            }
            catch ( final Exception e ) {
                throw new NHXFormatException( "malformed NHX format for event [" + nhx + "]:" + e.getMessage() );
            }
        }
    }

    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        if ( isUnassigned() ) {
        }
        else if ( isSpeciationOrDuplication() ) {
            sb.append( "?" );
        }
        else if ( isOther() || isRoot() || isTransfer() || isFusion() ) {
            sb.append( getEventType().toString() );
        }
        else {
            if ( getNumberOfDuplications() > 0 ) {
                if ( getNumberOfDuplications() > 1 ) {
                    sb.append( getNumberOfDuplications() );
                }
                sb.append( "D" );
            }
            if ( getNumberOfSpeciations() > 0 ) {
                if ( getNumberOfSpeciations() > 1 ) {
                    sb.append( getNumberOfSpeciations() );
                }
                sb.append( "S" );
            }
            if ( getNumberOfGeneLosses() > 0 ) {
                if ( getNumberOfGeneLosses() > 1 ) {
                    sb.append( getNumberOfGeneLosses() );
                }
                sb.append( "L" );
            }
        }
        return sb;
    }

    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        if ( isUnassigned() || isSpeciationOrDuplication() || isOther() || isRoot() || isTransfer() || isFusion() ) {
            sb.append( getEventType().toString() );
        }
        else {
            if ( isDuplication() ) {
                if ( getNumberOfDuplications() == 1 ) {
                    sb.append( "duplication" );
                }
                else {
                    sb.append( "duplications [" + getNumberOfDuplications() + "]" );
                }
            }
            else if ( isSpeciation() ) {
                if ( getNumberOfSpeciations() == 1 ) {
                    sb.append( "speciation" );
                }
                else {
                    sb.append( "speciations [" + getNumberOfSpeciations() + "]" );
                }
            }
            else if ( isGeneLoss() ) {
                if ( getNumberOfGeneLosses() == 1 ) {
                    sb.append( "gene-loss" );
                }
                else {
                    sb.append( "gene-losses [" + getNumberOfGeneLosses() + "]" );
                }
            }
            else {
                sb.append( "duplications [" + getNumberOfDuplications() + "] " );
                sb.append( "speciations [" + getNumberOfSpeciations() + "] " );
                sb.append( "gene-losses [" + getNumberOfGeneLosses() + "]" );
            }
        }
        return sb;
    }

    public PhylogenyData copy() {
        if ( isUnassigned() ) {
            return new Event();
        }
        else if ( _event_type != EventType.mixed ) {
            return new Event( _event_type );
        }
        else {
            return new Event( _duplications, _speciations, _gene_losses );
        }
    }

    public Confidence getConfidence() {
        return _confidence;
    }

    public EventType getEventType() {
        return _event_type;
    }

    public int getNumberOfDuplications() {
        return _duplications;
    }

    public int getNumberOfGeneLosses() {
        return _gene_losses;
    }

    public int getNumberOfSpeciations() {
        return _speciations;
    }

    /**
     * Returns true if this event contains one or more duplications events only
     * 
     * @return true if this event contains one or more duplications events only
     */
    public boolean isDuplication() {
        return ( _duplications > 0 ) && ( _gene_losses < 1 ) && ( _speciations < 1 );
    }

    public boolean isEqual( final PhylogenyData event ) {
        if ( ( event == null ) || !( event instanceof Event ) ) {
            return false;
        }
        final Event e = ( Event ) event;
        if ( getEventType().compareTo( e.getEventType() ) != 0 ) {
            return false;
        }
        if ( getNumberOfDuplications() != e.getNumberOfDuplications() ) {
            return false;
        }
        if ( getNumberOfSpeciations() != e.getNumberOfSpeciations() ) {
            return false;
        }
        if ( getNumberOfGeneLosses() != e.getNumberOfGeneLosses() ) {
            return false;
        }
        return true;
    }

    public boolean isFusion() {
        return _event_type == EventType.fusion;
    }

    /**
     * Returns true if this event contains one or more gene loss events only
     * 
     * @return true if this event contains one or more gene loss events only
     */
    public boolean isGeneLoss() {
        return ( _duplications < 1 ) && ( _gene_losses > 0 ) && ( _speciations < 1 );
    }

    public boolean isOther() {
        return _event_type == EventType.other;
    }

    public boolean isRoot() {
        return _event_type == EventType.root;
    }

    /**
     * Returns true if this event contains one or more speciation events only
     * 
     * @return true if this event contains one or more speciation events only
     */
    public boolean isSpeciation() {
        return ( _duplications < 1 ) && ( _gene_losses < 1 ) && ( _speciations > 0 );
    }

    public boolean isSpeciationOrDuplication() {
        return _event_type == EventType.speciation_or_duplication;
    }

    public boolean isTransfer() {
        return _event_type == EventType.transfer;
    }

    public boolean isUnassigned() {
        return ( _duplications == DEFAULT_VALUE ) && ( _event_type == EventType.unassigned );
    }

    public void setConfidence( final Confidence confidence ) {
        _confidence = confidence;
    }

    public void setDuplications( final int duplications ) {
        _duplications = duplications;
        _event_type = EventType.mixed;
    }

    public void setGeneLosses( final int gene_losses ) {
        _gene_losses = gene_losses;
        _event_type = EventType.mixed;
    }

    public void setSpeciations( final int speciations ) {
        _speciations = speciations;
        _event_type = EventType.mixed;
    }

    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        if ( !isUnassigned() && ( isSpeciationOrDuplication() || isDuplication() || isSpeciation() ) ) {
            sb.append( ":" );
            sb.append( NHXtags.IS_DUPLICATION );
            if ( isSpeciationOrDuplication() ) {
                sb.append( "?" );
            }
            else if ( isDuplication() ) {
                sb.append( "Y" );
            }
            else if ( isSpeciation() ) {
                sb.append( "N" );
            }
        }
        return sb;
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.EVENTS );
        if ( ( getEventType() != EventType.unassigned ) && ( getEventType() != EventType.mixed ) ) {
            PhylogenyDataUtil
                    .appendElement( writer, PhyloXmlMapping.EVENT_TYPE, getEventType().toString(), indentation );
        }
        if ( getNumberOfDuplications() > 0 ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.EVENT_DUPLICATIONS,
                                             getNumberOfDuplications() + "",
                                             indentation );
        }
        if ( getNumberOfSpeciations() > 0 ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.EVENT_SPECIATIONS,
                                             getNumberOfSpeciations() + "",
                                             indentation );
        }
        if ( getNumberOfGeneLosses() > 0 ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.EVENT_LOSSES,
                                             getNumberOfGeneLosses() + "",
                                             indentation );
        }
        if ( getConfidence() != null ) {
            getConfidence().toPhyloXML( writer, level, indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.EVENTS );
    }

    @Override
    public String toString() {
        return asText().toString();
    }

    public static Event createSingleDuplicationEvent() {
        return new Event( 1, 0, 0 );
    }

    public static Event createSingleSpeciationEvent() {
        return new Event( 0, 1, 0 );
    }

    public static Event createSingleSpeciationOrDuplicationEvent() {
        return new Event( EventType.speciation_or_duplication );
    }

    public static enum EventType {
        transfer, fusion, root, speciation_or_duplication, other, mixed, unassigned
    }
}
