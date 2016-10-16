/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <shark/future.hpp>
#include <shark/globals.hpp>

#include <cassert>

using namespace std;
using namespace shark;

std::set<std::shared_ptr<Handle>> Handle::ongoing;

void Handle::make_progress() {
    for (auto it = ongoing.begin(); it != ongoing.end(); ) {
        if ((*it)->test()) it = ongoing.erase(it);
        else ++it;
    }
}

Handle::Handle() {
}

Handle::~Handle() {
}

DoneHandle::DoneHandle() {
}

DoneHandle::~DoneHandle() {
}

bool DoneHandle::test() {
	return true;
}

void DoneHandle::wait() {
}

template<typename T>
Future<T>::Future(): done(false), h(), val() {
}

Future<void>::Future(): done(false), h() {
}

template<typename T>
Future<T>::~Future() {
	assert(!h || done);
}

Future<void>::~Future() {
	assert(!h || done);
}

template<typename T>
Future<T>::Future(Future<T>&& f): done(f.done), h(f.h), val(std::move(f.val)) {
}

Future<void>::Future(Future<void>&& f): done(f.done), h(f.h) {
}

template<typename T>
Future<T>& Future<T>::operator=(Future<T>&& f) {
	done = f.done;
        f.done = true;
	h = f.h;
	val = std::move(f.val);
	return *this;
}

Future<void>& Future<void>::operator=(Future<void>&& f) {
	done = f.done;
        f.done = true;
	h = f.h;
	return *this;
}

template<typename T>
Future<T>::Future(shared_ptr<Handle> h): done(false), h(h), val(make_unique<T>()) {
    Handle::ongoing.insert(h);
}

template<typename T>
Future<T>::Future(shared_ptr<Handle> h, unique_ptr<T>&& val): done(false), h(h), val(std::move(val)) {
    Handle::ongoing.insert(h);
}

Future<void>::Future(shared_ptr<Handle> h): done(false), h(h) {
    Handle::ongoing.insert(h);
}

template<typename T>
bool Future<T>::test() {
	return done || (done = h->test());
}

bool Future<void>::test() {
	return done || (done = h->test());
}

template<typename T>
const T& Future<T>::wait() {
	if(!done) {
		h->wait();
		done = true;
	}
	return *val;
}

void Future<void>::wait() {
	if(!done) {
		h->wait();
		done = true;
	}
}

// Set-up instantiations

#include "comm_impl.hpp"

#define SYMBT(T) template class shark::Future<T>;
#include "comm_int_inst"
#include "comm_fp_inst"
#include "comm_cplx_inst"
#include "comm_other_inst"
#undef SYMBT
