/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <shark/future.hpp>

#include <cassert>

using namespace std;
using namespace shark;

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
Future<T>::Future(Future<T>&& f): done(f.done), h(std::move(f.h)), val(std::move(f.val)) {
}

Future<void>::Future(Future<void>&& f): done(f.done), h(std::move(f.h)) {
}

template<typename T>
Future<T>& Future<T>::operator=(Future<T>&& f) {
	done = f.done;
	h = std::move(f.h);
	val = std::move(f.val);
	return *this;
}

Future<void>& Future<void>::operator=(Future<void>&& f) {
	done = f.done;
	h = std::move(f.h);
	return *this;
}

template<typename T>
Future<T>::Future(unique_ptr<Handle>&& h): done(false), h(std::move(h)), val(make_unique<T>()) {
}

template<typename T>
Future<T>::Future(unique_ptr<Handle>&& h, unique_ptr<T>&& val): done(false), h(std::move(h)), val(std::move(val)) {
}

Future<void>::Future(unique_ptr<Handle>&& h): done(false), h(std::move(h)) {
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
