// $Id: SimpleDomain.java,v 1.5 2009/11/17 03:51:34 cmzmasek Exp $
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

/*
 * A limited implementation of Domain. Its intended use is for when only a
 * domain identifier is needed. Note intended for general use.
 */
public class SimpleDomain implements Domain {

    final private DomainId _id;

    public SimpleDomain( final String id_str ) {
        if ( ForesterUtil.isEmpty( id_str ) ) {
            throw new IllegalArgumentException( "attempt to create protein domain with null or empty id" );
        }
        _id = new DomainId( id_str );
    }

    @Override
    public void addGoId( final GoId go_id ) {
        throw new RuntimeException( "method not implemented" );
    }

    public int compareTo( final Domain domain ) {
        if ( this == domain ) {
            return 0;
        }
        return getDomainId().compareTo( domain.getDomainId() );
    }

    public DomainId getDomainId() {
        return _id;
    }

    public int getFrom() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public GoId getGoId( final int i ) {
        throw new RuntimeException( "method not implemented" );
    }

    public int getLength() {
        throw new RuntimeException( "method not implemented" );
    }

    public short getNumber() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public int getNumberOfGoIds() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public double getPerDomainEvalue() {
        throw new RuntimeException( "method not implemented" );
    }

    @Override
    public double getPerDomainScore() {
        throw new RuntimeException( "method not implemented" );
    }

    public double getPerSequenceEvalue() {
        throw new RuntimeException( "method not implemented" );
    }

    public double getPerSequenceScore() {
        throw new RuntimeException( "method not implemented" );
    }

    public String getSearchParameter() {
        throw new RuntimeException( "method not implemented" );
    }

    public int getTo() {
        throw new RuntimeException( "method not implemented" );
    }

    public short getTotalCount() {
        throw new RuntimeException( "method not implemented" );
    }

    public boolean isCompleteQueryMatch() {
        throw new RuntimeException( "method not implemented" );
    }

    public boolean isCompleteTargetMatch() {
        throw new RuntimeException( "method not implemented" );
    }
}
