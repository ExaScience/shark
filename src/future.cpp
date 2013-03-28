/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
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
Future<T>::Future(): done(false), val() {
}

Future<void>::Future(): done(false), h() {
}

template<typename T>
Future<T>::~Future() {
	assert(!val || done);
}

Future<void>::~Future() {
	assert(!h || done);
}

template<typename T>
Future<T>::Future(Future<T>&& f): done(f.done), val(std::move(f.val)) {
}

Future<void>::Future(Future<void>&& f): done(f.done), h(std::move(f.h)) {
}

template<typename T>
Future<T>& Future<T>::operator=(Future<T>&& f) {
	done = f.done;
	val = std::move(f.val);
	return *this;
}

Future<void>& Future<void>::operator=(Future<void>&& f) {
	done = f.done;
	h = std::move(f.h);
	return *this;
}

template<typename T>
Future<T>::Future(bool done): done(done), val(new T()) {
}

Future<void>::Future(bool done, Handle* h): done(done), h(h) {
}

template<typename T>
Future<T>::Future(bool done, const T& val): done(done), val(new T(val)) {
}

template<typename T>
Future<T>::Future(bool done, T&& val): done(done), val(new T(std::move(val))) {
}

bool Future<void>::test() {
	return done || done = h->test();
}

void Future<void>::wait() {
	if(done)
		return;
	h->wait();
	done = true;
}

// Set-up instantiations

#include "types"

//template class Future<int>;
