// $Id: BasicDomain.java,v 1.9 2009/11/17 03:51:34 cmzmasek Exp $
//
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

package org.forester.surfacing;

import org.forester.go.GoId;
import org.forester.util.ForesterUtil;

public class BasicDomain implements Domain {

    final private DomainId _id;
    final private int      _from;
    final private int      _to;
    final private short    _number;
    final private short    _total_count;
    final private double   _per_sequence_evalue;
    final private double   _per_sequence_score;
    final private double   _per_domain_evalue;
    final private double   _per_domain_score;

    public BasicDomain( final String id_str ) {
        if ( ForesterUtil.isEmpty( id_str ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain with null or empty id" );
        }
        _id = new DomainId( id_str );
        _from = -1;
        _to = -1;
        _number = -1;
        _total_count = -1;
        _per_sequence_evalue = -1;
        _per_sequence_score = -1;
        _per_domain_evalue = -1;
        _per_domain_score = -1;
    }

    public BasicDomain( final String id_str,
                        final int from,
                        final int to,
                        final short number,
                        final short total_count,
                        final double per_sequence_evalue,
                        final double per_sequence_score ) {
        if ( ( from >= to ) || ( from < 0 ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain from " + from + " to " + to );
        }
        if ( ForesterUtil.isEmpty( id_str ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain with null or empty id" );
        }
        if ( ( number > total_count ) || ( number < 0 ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain number " + number + " out of "
                    + total_count );
        }
        if ( per_sequence_evalue < 0.0 ) {
            throw new IllegalArgumentException( "attempt to create protein domain with E-value" );
        }
        _id = new DomainId( id_str );
        _from = from;
        _to = to;
        _number = number;
        _total_count = total_count;
        _per_sequence_evalue = per_sequence_evalue;
        _per_sequence_score = per_sequence_score;
        _per_domain_evalue = -1;
        _per_domain_score = -1;
    }

    public BasicDomain( final String id_str,
                        final int from,
                        final int to,
                        final short number,
                        final short total_count,
                        final double per_sequence_evalue,
                        final double per_sequence_score,
                        final double per_domain_evalue,
                        final double per_domain_score ) {
        if ( ( from >= to ) || ( from < 0 ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain from " + from + " to " + to );
        }
        if ( ForesterUtil.isEmpty( id_str ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain with null or empty id" );
        }
        if ( ( number > total_count ) || ( number < 0 ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain number " + number + " out of "
                    + total_count );
        }
        if ( ( per_sequence_evalue < 0.0 ) || ( per_domain_evalue < 0.0 ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain with E-value" );
        }
        _id = new DomainId( id_str );
        _from = from;
        _to = to;
        _number = number;
        _total_count = total_count;
        _per_sequence_evalue = per_sequence_evalue;
        _per_sequence_score = per_sequence_score;
        _per_domain_evalue = per_domain_evalue;
        _per_domain_score = per_domain_score;
    }

    public void addGoId( final GoId go_id ) {
        getDomainId().getGoIds().add( go_id );
    }

    /**
     * Basic domains are compared/sorted based upon their identifiers (case
     * insensitive) and their numbers.
     * 
     */
    public int compareTo( final Domain domain ) {
        if ( domain.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to compare [" + domain.getClass() + "] to " + "["
                    + this.getClass() + "]" );
        }
        if ( this == domain ) {
            return 0;
        }
        return getDomainId().compareTo( domain.getDomainId() );
    }

    /**
     * Basic domains are considered equal if they have the same identifier (case
     * sensitive).
     * 
     */
    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check [" + this.getClass() + "] equality to " + o + " ["
                    + o.getClass() + "]" );
        }
        else {
            return getDomainId().equals( ( ( Domain ) o ).getDomainId() );
        }
    }

    public DomainId getDomainId() {
        return _id;
    }

    public int getFrom() {
        return _from;
    }

    public GoId getGoId( final int i ) {
        return getDomainId().getGoIds().get( i );
    }

    public short getNumber() {
        return _number;
    }

    public int getNumberOfGoIds() {
        return getDomainId().getGoIds().size();
    }

    @Override
    public double getPerDomainEvalue() {
        return _per_domain_evalue;
    }

    @Override
    public double getPerDomainScore() {
        return _per_domain_score;
    }

    public double getPerSequenceEvalue() {
        return _per_sequence_evalue;
    }

    public double getPerSequenceScore() {
        return _per_sequence_score;
    }

    public int getTo() {
        return _to;
    }

    public short getTotalCount() {
        return _total_count;
    }

    @Override
    public int hashCode() {
        return getDomainId().getId().hashCode();
    }

    @Override
    public String toString() {
        return toStringBuffer().toString();
    }

    public StringBuffer toStringBuffer() {
        return new StringBuffer( getDomainId().getId() );
    }
}
