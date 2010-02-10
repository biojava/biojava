// $Id: Sequence.java,v 1.57 2009/12/12 00:14:39 cmzmasek Exp $
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
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlDataFormatException;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.parsers.phyloxml.PhyloXmlUtil;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class Sequence implements PhylogenyData {

    private String              _mol_sequence;
    private String              _name;
    private Accession           _accession;
    private String              _symbol;
    private String              _location;
    private String              _type;
    private List<PhylogenyData> _annotations;
    private DomainArchitecture  _da;
    private Uri                 _uri;

    public Sequence() {
        init();
    }

    public void addAnnotation( final Annotation annotation ) {
        getAnnotations().add( annotation );
    }

    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        if ( getAccession() != null ) {
            sb.append( "[" );
            sb.append( getAccession() );
            sb.append( "] " );
        }
        if ( !ForesterUtil.isEmpty( getName() ) ) {
            sb.append( getName() );
            sb.append( " " );
        }
        if ( !ForesterUtil.isEmpty( getLocation() ) ) {
            sb.append( getLocation() );
        }
        return sb;
    }

    public StringBuffer asText() {
        return asSimpleText();
    }

    /**
     * Not a deep copy.
     * 
     */
    public PhylogenyData copy() {
        final Sequence seq = new Sequence();
        seq.setAnnotations( getAnnotations() );
        seq.setName( new String( getName() ) );
        seq.setSymbol( new String( getSymbol() ) );
        seq.setMolecularSequence( new String( getMolecularSequence() ) );
        seq.setLocation( new String( getLocation() ) );
        if ( getAccession() != null ) {
            seq.setAccession( ( Accession ) getAccession().copy() );
        }
        else {
            seq.setAccession( null );
        }
        seq.setType( new String( getType() ) );
        if ( getUri() != null ) {
            seq.setUri( ( Uri ) getUri().copy() );
        }
        else {
            seq.setUri( null );
        }
        if ( getDomainArchitecture() != null ) {
            seq.setDomainArchitecture( ( DomainArchitecture ) getDomainArchitecture().copy() );
        }
        else {
            seq.setDomainArchitecture( null );
        }
        return seq;
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            return false;
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return isEqual( ( Sequence ) o );
        }
    }

    public Accession getAccession() {
        return _accession;
    }

    public PhylogenyData getAnnotation( final int index ) {
        return getAnnotations().get( index );
    }

    public List<PhylogenyData> getAnnotations() {
        return _annotations;
    }

    public DomainArchitecture getDomainArchitecture() {
        return _da;
    }

    public String getLocation() {
        return _location;
    }

    public String getMolecularSequence() {
        return _mol_sequence;
    }

    public String getName() {
        return _name;
    }

    public String getSymbol() {
        return _symbol;
    }

    public String getType() {
        return _type;
    }

    public Uri getUri() {
        return _uri;
    }

    @Override
    public int hashCode() {
        if ( getAccession() != null ) {
            return getAccession().hashCode();
        }
        int result = getSymbol().hashCode();
        if ( getName().length() > 0 ) {
            result ^= getName().hashCode();
        }
        if ( getMolecularSequence().length() > 0 ) {
            result ^= getMolecularSequence().hashCode();
        }
        return result;
    }

    public void init() {
        setAnnotations( new ArrayList<PhylogenyData>() );
        setName( "" );
        setMolecularSequence( "" );
        setLocation( "" );
        setAccession( null );
        setSymbol( "" );
        setType( "" );
        setDomainArchitecture( null );
        setUri( null );
    }

    public boolean isEmpty() {
        return ( getAccession() == null ) && ForesterUtil.isEmpty( getName() ) && ForesterUtil.isEmpty( getSymbol() )
                && ForesterUtil.isEmpty( getType() ) && ForesterUtil.isEmpty( getLocation() )
                && ForesterUtil.isEmpty( getMolecularSequence() ) && ( getDomainArchitecture() == null )
                && getAnnotations().isEmpty() && ( getUri() == null );
    }

    public boolean isEqual( final PhylogenyData data ) {
        if ( this == data ) {
            return true;
        }
        final Sequence s = ( Sequence ) data;
        if ( ( getAccession() != null ) && ( s.getAccession() != null ) ) {
            return getAccession().isEqual( s.getAccession() );
        }
        return s.getMolecularSequence().equals( getMolecularSequence() ) && s.getName().equals( getName() )
                && s.getSymbol().equals( getSymbol() );
    }

    public void setAccession( final Accession accession ) {
        _accession = accession;
    }

    private void setAnnotations( final List<PhylogenyData> annotations ) {
        _annotations = annotations;
    }

    public void setDomainArchitecture( final DomainArchitecture ds ) {
        _da = ds;
    }

    public void setLocation( final String description ) {
        _location = description;
    }

    public void setMolecularSequence( final String mol_sequence ) {
        _mol_sequence = mol_sequence;
    }

    public void setName( final String name ) {
        _name = name;
    }

    public void setSymbol( final String symbol ) {
        if ( !ForesterUtil.isEmpty( symbol ) && !PhyloXmlUtil.SEQUENCE_SYMBOL_PATTERN.matcher( symbol ).matches() ) {
            throw new PhyloXmlDataFormatException( "illegal sequence symbol: [" + symbol + "]" );
        }
        _symbol = symbol;
    }

    public void setType( final String type ) {
        if ( !ForesterUtil.isEmpty( type ) && !PhyloXmlUtil.SEQUENCE_TYPES.contains( type ) ) {
            throw new PhyloXmlDataFormatException( "illegal sequence type: [" + type + "]" );
        }
        _type = type;
    }

    public void setUri( final Uri uri ) {
        _uri = uri;
    }

    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        if ( getName().length() > 0 ) {
            sb.append( ":" );
            sb.append( NHXtags.GENE_NAME );
            sb.append( ForesterUtil.replaceIllegalNhxCharacters( getName() ) );
        }
        if ( getAccession() != null ) {
            getAccession().toNHX();
        }
        if ( getDomainArchitecture() != null ) {
            sb.append( getDomainArchitecture().toNHX() );
        }
        return sb;
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        final String my_ind = indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE;
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer, PhyloXmlMapping.SEQUENCE, PhyloXmlMapping.SEQUENCE_TYPE, getType() );
        if ( !ForesterUtil.isEmpty( getSymbol() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.SEQUENCE_SYMBOL, getSymbol(), indentation );
        }
        if ( ( getAccession() != null ) && !ForesterUtil.isEmpty( getAccession().getValue() ) ) {
            getAccession().toPhyloXML( writer, level, indentation );
        }
        if ( !ForesterUtil.isEmpty( getName() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.SEQUENCE_NAME, getName(), indentation );
        }
        if ( !ForesterUtil.isEmpty( getLocation() ) ) {
            PhylogenyDataUtil.appendElement( writer, PhyloXmlMapping.SEQUENCE_LOCATION, getLocation(), indentation );
        }
        if ( !ForesterUtil.isEmpty( getMolecularSequence() ) ) {
            PhylogenyDataUtil.appendElement( writer,
                                             PhyloXmlMapping.SEQUENCE_MOL_SEQ,
                                             getMolecularSequence(),
                                             indentation );
        }
        if ( getUri() != null ) {
            getUri().toPhyloXML( writer, level, indentation );
        }
        if ( !getAnnotations().isEmpty() ) {
            for( final PhylogenyData annotation : getAnnotations() ) {
                annotation.toPhyloXML( writer, level, my_ind );
            }
        }
        if ( getDomainArchitecture() != null ) {
            getDomainArchitecture().toPhyloXML( writer, level, my_ind );
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.SEQUENCE );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
