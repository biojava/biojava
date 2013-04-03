// $Id: DomainArchitecture.java,v 1.17 2009/10/26 23:29:39 cmzmasek Exp $
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
import java.util.SortedMap;
import java.util.StringTokenizer;
import java.util.TreeMap;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.io.writers.PhylogenyWriter;
import org.forester.util.ForesterUtil;

public class DomainArchitecture implements PhylogenyData {

    public final static String               NHX_SEPARATOR = ">";
    private static final double              INCREASE_KEY  = 0.0001;
    private SortedMap<Double, ProteinDomain> _domains;
    private int                              _total_length;

    public DomainArchitecture() {
        init();
    }

    public DomainArchitecture( final List<PhylogenyData> domains, final int total_length ) {
        init();
        for( final PhylogenyData phylogenyData : domains ) {
            final ProteinDomain pd = ( ProteinDomain ) phylogenyData;
            addDomain( pd );
        }
        _total_length = total_length;
    }

    public DomainArchitecture( final String da_str ) {
        init();
        int total_length = 0;
        int to = -1;
        try {
            final StringTokenizer st = new StringTokenizer( da_str, DomainArchitecture.NHX_SEPARATOR );
            final String length_str = ( String ) st.nextElement();
            total_length = new Integer( length_str ).intValue();
            while ( st.hasMoreElements() ) {
                final String from_str = ( String ) st.nextElement();
                final String to_str = ( String ) st.nextElement();
                final String support_str = ( String ) st.nextElement();
                final String name = ( String ) st.nextElement();
                to = new Integer( to_str ).intValue();
                final int from = new Integer( from_str ).intValue();
                final double support = new Double( support_str ).doubleValue();
                final ProteinDomain pd = new ProteinDomain( name, from, to, support );
                addDomain( pd );
            }
        }
        catch ( final Exception e ) {
            throw new IllegalArgumentException( "Malformed format for domain structure \"" + da_str + "\": "
                    + e.getMessage() );
        }
        if ( to > total_length ) {
            throw new IllegalArgumentException( "total length of domain structure is too short" );
        }
        _total_length = total_length;
    }

    public void addDomain( final ProteinDomain pd ) {
        Double key = new Double( pd.getFrom() );
        while ( _domains.containsKey( key ) ) {
            key = new Double( key.doubleValue() + DomainArchitecture.INCREASE_KEY );
        }
        _domains.put( key, pd );
    }

    public StringBuffer asSimpleText() {
        final StringBuffer sb = new StringBuffer();
        for( int i = 0; i < getDomains().size(); ++i ) {
            if ( i > 0 ) {
                sb.append( "~" );
            }
            sb.append( getDomain( i ).asSimpleText() );
        }
        return sb;
    }

    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        for( int i = 0; i < getDomains().size(); ++i ) {
            if ( i > 0 ) {
                sb.append( "~" );
            }
            sb.append( getDomain( i ).asText() );
        }
        return sb;
    }

    public PhylogenyData copy() {
        final List<PhylogenyData> domains = new ArrayList<PhylogenyData>( getDomains().size() );
        for( int i = 0; i < getDomains().size(); ++i ) {
            domains.add( getDomain( i ).copy() );
        }
        return new DomainArchitecture( domains, getTotalLength() );
    }

    public ProteinDomain getDomain( final int i ) {
        return ( ProteinDomain ) _domains.values().toArray()[ i ];
    }

    public SortedMap<Double, ProteinDomain> getDomains() {
        return _domains;
    }

    public int getNumberOfDomains() {
        return _domains.size();
    }

    public int getTotalLength() {
        return _total_length;
    }

    private void init() {
        _domains = new TreeMap<Double, ProteinDomain>();
        _total_length = 0;
    }

    /**
     * Returns true if the names and the order of the domains match (domain and
     * linker lengths are ignored).
     * 
     * 
     */
    public boolean isEqual( final PhylogenyData domain_architecture ) {
        if ( domain_architecture == null ) {
            return false;
        }
        if ( !( domain_architecture instanceof DomainArchitecture ) ) {
            return false;
        }
        final DomainArchitecture d = ( DomainArchitecture ) domain_architecture;
        if ( getDomains().size() != d.getDomains().size() ) {
            return false;
        }
        for( int i = 0; i < getDomains().size(); ++i ) {
            if ( !getDomain( i ).getName().equals( d.getDomain( i ).getName() ) ) {
                return false;
            }
        }
        return true;
    }

    public void setTotalLength( final int total_length ) {
        _total_length = total_length;
    }

    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        sb.append( ":" );
        sb.append( NHXtags.DOMAIN_STRUCTURE );
        sb.append( getTotalLength() );
        if ( getDomains() != null ) {
            for( int i = 0; i < getDomains().size(); ++i ) {
                sb.append( DomainArchitecture.NHX_SEPARATOR );
                sb.append( getDomain( i ).getFrom() );
                sb.append( DomainArchitecture.NHX_SEPARATOR );
                sb.append( getDomain( i ).getTo() );
                sb.append( DomainArchitecture.NHX_SEPARATOR );
                sb.append( getDomain( i ).getConfidence() );
                sb.append( DomainArchitecture.NHX_SEPARATOR );
                sb.append( ForesterUtil.replaceIllegalNhxCharacters( getDomain( i ).getName() ) );
            }
        }
        return sb;
    }

    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendOpen( writer,
                                      PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECURE,
                                      PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECTURE_LENGTH,
                                      getTotalLength() + "" );
        if ( getDomains() != null ) {
            for( int i = 0; i < getDomains().size(); ++i ) {
                getDomain( i ).toPhyloXML( writer, level, indentation + PhylogenyWriter.PHYLO_XML_INTENDATION_BASE );
            }
        }
        writer.write( ForesterUtil.LINE_SEPARATOR );
        writer.write( indentation );
        PhylogenyDataUtil.appendClose( writer, PhyloXmlMapping.SEQUENCE_DOMAIN_ARCHITECURE );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
