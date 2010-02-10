// $Id: MappingResults.java,v 1.1 2008/06/26 22:48:46 cmzmasek Exp $
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
// WWW: www.phylosoft.org

package org.forester.surfacing;

public class MappingResults {

    private String _description;
    private int    _sum_of_successes;
    private int    _sum_of_failures;

    public String getDescription() {
        return _description;
    }

    public int getSumOfFailures() {
        return _sum_of_failures;
    }

    public int getSumOfSuccesses() {
        return _sum_of_successes;
    }

    public void setDescription( final String description ) {
        _description = description;
    }

    public void setSumOfFailures( final int sum_of_failures ) {
        _sum_of_failures = sum_of_failures;
    }

    public void setSumOfSuccesses( final int sum_of_successes ) {
        _sum_of_successes = sum_of_successes;
    }
}
