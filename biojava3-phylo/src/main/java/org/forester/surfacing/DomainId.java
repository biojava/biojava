// $Id: DomainId.java,v 1.8 2009/10/26 23:29:40 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.List;

import org.forester.go.GoId;
import org.forester.util.ForesterUtil;

public class DomainId implements Comparable<DomainId> {

    final private String _id;
    private List<GoId>   _go_ids;

    public DomainId( final String id ) {
        if ( ForesterUtil.isEmpty( id ) ) {
            throw new IllegalArgumentException( "attempt to create domain id from empty or null string" );
        }
        _id = id.trim();
        if ( _id.indexOf( ' ' ) > -1 ) {
            throw new IllegalArgumentException( "attempt to create domain id from string containing one ore more spaces ["
                    + _id + "]" );
        }
        else if ( _id.indexOf( BinaryDomainCombination.SEPARATOR ) > -1 ) {
            throw new IllegalArgumentException( "attempt to create domain id from string containing the separator character ["
                    + BinaryDomainCombination.SEPARATOR + "] for domain combinations [" + _id + "]" );
        }
        setGoIds( null );
    }

    public void addGoId( final GoId go_id ) {
        if ( getGoIds() == null ) {
            setGoIds( new ArrayList<GoId>() );
        }
        getGoIds().add( go_id );
    }

    @Override
    public int compareTo( final DomainId domain_id ) {
        if ( this == domain_id ) {
            return 0;
        }
        return getId().toLowerCase().compareTo( domain_id.getId().toLowerCase() );
    }

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
            return getId().equals( ( ( DomainId ) o ).getId() );
        }
    }

    public GoId getGoId( final int i ) {
        return getGoIds().get( i );
    }

    // Note.
    // The fact that equals and compareTo do not behave the same in cases where ids only differ by their case
    // is not ideal. From Sun regarding Interface SortedSet<E>:
    // "Note that the ordering maintained by a sorted set (whether or not an explicit comparator is provided)
    // must be consistent with equals if the sorted set is to correctly implement the Set interface.
    // (See the Comparable interface or Comparator interface for a precise definition of consistent 
    // with equals.) This is so because the Set interface is defined in terms of the equals  operation,
    // but a sorted set performs all element comparisons using its compareTo (or compare) method, 
    // so two elements that are deemed equal by this method are, from the standpoint of the sorted set,
    // equal. The behavior of a sorted set is well-defined even if its ordering is inconsistent with equals; 
    // it just fails to obey the general contract of the Set interface."
    List<GoId> getGoIds() {
        return _go_ids;
    }

    public String getId() {
        return _id;
    }

    public int getNumberOfGoIds() {
        if ( getGoIds() == null ) {
            return 0;
        }
        return getGoIds().size();
    }

    @Override
    public int hashCode() {
        return getId().hashCode();
    }

    private void setGoIds( final List<GoId> go_ids ) {
        _go_ids = go_ids;
    }

    @Override
    public String toString() {
        return getId();
    }
}
